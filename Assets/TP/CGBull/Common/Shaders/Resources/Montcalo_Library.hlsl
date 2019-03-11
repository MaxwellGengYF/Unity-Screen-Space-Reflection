#ifndef _Montcalo_Library_
#define _Montcalo_Library_

#include "Common.hlsl"

uint ReverseBits32(uint bits)
{
	bits = (bits << 16) | (bits >> 16);
	bits = ((bits & 0x00ff00ff) << 8) | ((bits & 0xff00ff00) >> 8);
	bits = ((bits & 0x0f0f0f0f) << 4) | ((bits & 0xf0f0f0f0) >> 4);
	bits = ((bits & 0x33333333) << 2) | ((bits & 0xcccccccc) >> 2);
	bits = ((bits & 0x55555555) << 1) | ((bits & 0xaaaaaaaa) >> 1);
	return bits;
}

uint2 SobolIndex(uint2 Base, int Index, int Bits = 10) {
	uint2 SobolNumbers[10] = {
		uint2(0x8680u, 0x4c80u), uint2(0xf240u, 0x9240u), uint2(0x8220u, 0x0e20u), uint2(0x4110u, 0x1610u), uint2(0xa608u, 0x7608u),
		uint2(0x8a02u, 0x280au), uint2(0xe204u, 0x9e04u), uint2(0xa400u, 0x4682u), uint2(0xe300u, 0xa74du), uint2(0xb700u, 0x9817u),
	};

	uint2 Result = Base;
	[ROLL] 
    for (int b = 0; b < 10 && b < Bits; ++b) {
		Result ^= (Index & (1 << b)) ? SobolNumbers[b] : 0;
	}
	return Result;
}

uint HaltonSequence(uint Index, uint base = 3)
{
	uint result = 0;
	uint f = 1;
	uint i = Index;
	
	UNITY_UNROLL
	while (i > 0) {
		f = f / base;
		result = result + f * (i % base);
		i = floor(i / base);
	}
	return result;
}

float2 Hammersley(uint Index, uint NumSamples)
{
	return float2((float)Index / (float)NumSamples, ReverseBits32(Index));
}

float2 Hammersley(uint Index, uint NumSamples, uint2 Random)
{
	float E1 = frac((float)Index / NumSamples + float(Random.x & 0xffff) / (1 << 16));
	float E2 = float(ReverseBits32(Index) ^ Random.y) * 2.3283064365386963e-10;
	return float2(E1, E2);
}

float3x3 GetTangentBasis(float3 TangentZ) {
	float3 UpVector = abs(TangentZ.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
	float3 TangentX = normalize(cross( UpVector, TangentZ));
	float3 TangentY = cross(TangentZ, TangentX);
	return float3x3(TangentX, TangentY, TangentZ);
}

float3 TangentToWorld(float3 Vec, float3 TangentZ)
{
	return mul(Vec, GetTangentBasis(TangentZ));
}

float4 TangentToWorld(float3 Vec, float4 TangentZ)
{
	half3 T2W = TangentToWorld(Vec, TangentZ.rgb);
	return half4(T2W, TangentZ.a);
}

float2 RandToCircle(uint2 Rand) {
	float2 sf = float2(Rand) * (sqrt(2.) / 0xffff) - sqrt(0.5);	
	float2 sq = sf*sf;
	float root = sqrt(2.*max(sq.x, sq.y) - min(sq.x, sq.y));
	if (sq.x > sq.y) {
		sf.x = sf.x > 0 ? root : -root;
	}
	else {
		sf.y = sf.y > 0 ? root : -root;
	}
	return sf;
}

float4 UniformSampleSphere(float2 E) {
	float Phi = 2 * PI * E.x;
	float CosTheta = 1 - 2 * E.y;
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;

	float PDF = 1 / (4 * PI);

	return float4(H, PDF);
}

float4 UniformSampleHemisphere(float2 E) {
	float Phi = 2 * PI * E.x;
	float CosTheta = E.y;
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 H;
	H.x = SinTheta * cos( Phi );
	H.y = SinTheta * sin( Phi );
	H.z = CosTheta;

	float PDF = 1.0 / (2 * PI);
	return float4(H, PDF);
}

float2 UniformSampleDisk(float2 Random) {
	const float Theta = 2.0f * (float)PI * Random.x;
	const float Radius = sqrt(Random.y);
	return float2(Radius * cos(Theta), Radius * sin(Theta));
}

float4 CosineSampleHemisphere(float2 E) {
	float Phi = 2 * PI * E.x;
	float CosTheta = sqrt(E.y);
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;

	float PDF = CosTheta / PI;
	return float4(H, PDF);
}

float4 CosineSampleHemisphere(float2 E, float3 Normal) {
	float4 Sampler = CosineSampleHemisphere(E);
	return float4(normalize(Normal + Sampler.rgb), Sampler.a);
}

float4 UniformSampleCone(float2 E, float CosThetaMax) {
	float Phi = 2 * PI * E.x;
	float CosTheta = lerp(CosThetaMax, 1, E.y);
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 L;
	L.x = SinTheta * cos( Phi );
	L.y = SinTheta * sin( Phi );
	L.z = CosTheta;

	float PDF = 1.0 / (2 * PI * (1 - CosThetaMax));
	return float4(L, PDF);
}


float4 ImportanceSampleLambert(float2 E)
{
    float3 L = CosineSampleHemisphere(E).rgb;
	return float4(L, 1);
}

float4 ImportanceSampleBlinn(float2 E, float Roughness) {
	float m = Roughness * Roughness;
	float m2 = m * m;
		
	float Phi = 2 * PI * E.x;
	float n = 2 / m2 - 2;
	float CosTheta = pow(max(E.y, 0.001), 1 / (n + 1));
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;
		
	float D = (n + 2)/ (2 * PI) * saturate(pow(CosTheta, n));
	float pdf = D * CosTheta;
	return float4(H, pdf); 
}

float3 ImportanceSampleGGX(float2 E, float3 N, float Roughness)
{
	float a = Roughness * Roughness;

	float phi = 2.0 * PI * E.x;
	float cosTheta = sqrt((1.0 - E.y) / (1.0 + (a*a - 1.0) * E.y));
	float sinTheta = sqrt(1.0 - cosTheta*cosTheta);
	
	// from spherical coordinates to cartesian coordinates - halfway vector
	float3 H;
	H.x = cos(phi) * sinTheta;
	H.y = sin(phi) * sinTheta;
	H.z = cosTheta;
	
	// from tangent-space H vector to world-space sample vector
	float3 up          = abs(N.z) < 0.999 ? float3(0.0, 0.0, 1.0) : float3(1.0, 0.0, 0.0);
	float3 tangent   = normalize(cross(up, N));
	float3 bitangent = cross(N, tangent);
	
	float3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
	return normalize(sampleVec);
}

float4 ImportanceSampleGGX(float2 E, float Roughness) {
	float m = Roughness * Roughness;
	float m2 = m * m;

	float Phi = 2 * PI * E.x;
	float CosTheta = sqrt((1 - E.y) / ( 1 + (m2 - 1) * E.y));
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;
			
	float d = (CosTheta * m2 - CosTheta) * CosTheta + 1;
	float D = m2 / (PI * d * d);
			
	float PDF = D * CosTheta;

	return float4(H, PDF);
}

float4 ImportanceSampleInverseGGX(float2 E, float Roughness) {
	float m = Roughness * Roughness;
	float m2 = m * m;
	float A = 4;

	float Phi = 2 * PI * E.x;
	float CosTheta = sqrt((1 - E.y) / ( 1 + (m2 - 1) * E.y));
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	float3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;
			
	float d = (CosTheta - m2 * CosTheta) * CosTheta + m2;
	float D = rcp(Inv_PI * (1 + A * m2)) * (1 + 4 * m2 * m2 / (d * d));
			
	float PDF = D * CosTheta;

	return float4(H, PDF);
}

void SampleAnisoGGXDir(float2 u, float3 V, float3 N, float3 tX, float3 tY, float roughnessT, float roughnessB, out float3 H, out float3 L) {
    H = sqrt(u.x / (1 - u.x)) * (roughnessT * cos(Two_PI * u.y) * tX + roughnessB * sin(Two_PI * u.y) * tY) + N;
    H = normalize(H);
    L = 2 * saturate(dot(V, H)) * H - V;
}

void ImportanceSampleAnisoGGX(float2 u, float3 V, float3 N, float3 tX, float3 tY, float roughnessT, float roughnessB, float NoV, out float3 L, out float VoH, out float NoL, out float weightOverPdf)
{
    float3 H;
    SampleAnisoGGXDir(u, V, N, tX, tY, roughnessT, roughnessB, H, L);

    float NoH = saturate(dot(N, H));
    VoH = saturate(dot(V, H));
    NoL = saturate(dot(N, L));

    float ToV = dot(tX, V);
    float BoV = dot(tY, V);
    float ToL = dot(tX, L);
    float BoL = dot(tY, L);

    float aT = roughnessT;
    float aT2 = aT * aT;
    float aB = roughnessB;
    float aB2 = aB * aB;
    float lambdaV = NoL * sqrt(aT2 * ToV * ToV + aB2 * BoV * BoV + NoV * NoV);
    float lambdaL = NoV * sqrt(aT2 * ToL * ToL + aB2 * BoL * BoL + NoL * NoL);
    float Vis = 0.5 / (lambdaV + lambdaL);
	
    weightOverPdf = 4 * Vis * NoL * VoH / NoH;
}

float MISWeight(uint Num, float PDF, uint OtherNum, float OtherPDF) {
	float Weight = Num * PDF;
	float OtherWeight = OtherNum * OtherPDF;
	return Weight * Weight / (Weight * Weight + OtherWeight * OtherWeight);
}

#endif