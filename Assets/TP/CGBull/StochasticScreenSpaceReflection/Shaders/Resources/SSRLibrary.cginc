#include "UnityStandardBRDF.cginc"
#include "Assets/TP/CGBull/Common/Shaders/Resources/Include_HLSL.hlsl"

inline half GetScreenFadeBord(half2 pos, half value)
{
    half borderDist = min(1 - max(pos.x, pos.y), min(pos.x, pos.y));
    return saturate(borderDist > value ? 1 : borderDist / value);
}


/////////////////////////////////////LInear3DTrace/////////////////////////////////////
float4 LinearTraceRay3DSpace(Texture2D _DepthTexture, SamplerState sampler_DepthTexture, int NumSteps, float2 BlueNoise, float3 rayPos, float3 rayDir) {
	float mask = 0.0, endDepth = 0.0;
    float rayDepth = rayPos.z;

	float2 jitter = BlueNoise + 0.5;
	float StepSize = 1 / (float)NumSteps;
	StepSize = StepSize * (jitter.x + jitter.y) + StepSize;

    UNITY_LOOP
	for (int i = 0;  i < NumSteps; i++) {
		endDepth = Texture2DSampleLevel(_DepthTexture, sampler_DepthTexture, rayPos.xy, 0.0);
		if (rayDepth > endDepth) {
			mask = 1;
			break;
		}
		rayPos += rayDir * StepSize;
        rayDepth = rayPos.z;
	}
	return float4(rayPos, mask);
}


/////////////////////////////////////Linear2DTrace/////////////////////////////////////
inline half distanceSquared(half2 A, half2 B)
{
    A -= B;
    return dot(A, A);
}

inline half distanceSquared(half3 A, half3 B)
{
    A -= B;
    return dot(A, A);
}

void swap(inout half v0, inout half v1)
{
    half temp = v0;
    v0 = v1;
    v1 = temp;
}

bool intersectsDepthBuffer(half rayZMin, half rayZMax, half sceneZ, half layerThickness)
{
    return (rayZMax >= sceneZ - layerThickness) && (rayZMin <= sceneZ);
}

void rayIterations(sampler2D forntDepth, in bool traceBehind_Old, in bool traceBehind, inout half2 P, inout half stepDirection, inout half end, inout int stepCount, inout int maxSteps, inout bool intersecting,
                   inout half sceneZ, inout half2 dP, inout half3 Q, inout half3 dQ, inout half k, inout half dk,
                   inout half rayZMin, inout half rayZMax, inout half prevZMaxEstimate, inout bool permute, inout half2 hitPixel,
                   half2 invSize, inout half layerThickness)
{
    bool stop = intersecting;
    
    for (; (P.x * stepDirection) <= end && stepCount < maxSteps && !stop; P += dP, Q.z += dQ.z, k += dk, stepCount += 1)
    {
        rayZMin = prevZMaxEstimate;
        rayZMax = (dQ.z * 0.5 + Q.z) / (dk * 0.5 + k);
        prevZMaxEstimate = rayZMax;

        if (rayZMin > rayZMax) {
            swap(rayZMin, rayZMax);
        }

        hitPixel = permute ? P.yx : P;
        sceneZ = tex2Dlod(forntDepth, half4(hitPixel * invSize, 0, 0)).r;
        sceneZ = -LinearEyeDepth(sceneZ);
        bool isBehind = (rayZMin <= sceneZ);

        if (traceBehind_Old == 1) {
            intersecting = isBehind && (rayZMax >= sceneZ - layerThickness);
        } else {
            intersecting = (rayZMax >= sceneZ - layerThickness);
        }

        stop = traceBehind ? intersecting : isBehind;
    }
    P -= dP, Q.z -= dQ.z, k -= dk;
}

bool Linear2D_Trace(sampler2D forntDepth,
                             half3 csOrigin,
                             half3 csDirection,
                             half4x4 projectMatrix,
                             half2 csZBufferSize,
                             half jitter,
                             int maxSteps,
                             half layerThickness,
                             half traceDistance,
                             in out half2 hitPixel,
                             int stepSize,
                             bool traceBehind,
                             in out half3 csHitPoint,
                             in out half stepCount)
{

    half2 invSize = half2(1 / csZBufferSize.x, 1 / csZBufferSize.y);
    hitPixel = half2(-1, -1);

    half nearPlaneZ = -0.01;
    half rayLength = ((csOrigin.z + csDirection.z * traceDistance) > nearPlaneZ) ? ((nearPlaneZ - csOrigin.z) / csDirection.z) : traceDistance;
    half3 csEndPoint = csDirection * rayLength + csOrigin;
    half4 H0 = mul(projectMatrix, half4(csOrigin, 1));
    half4 H1 = mul(projectMatrix, half4(csEndPoint, 1));
    half k0 = 1 / H0.w;
    half k1 = 1 / H1.w;
    half2 P0 = H0.xy * k0;
    half2 P1 = H1.xy * k1;
    half3 Q0 = csOrigin * k0;
    half3 Q1 = csEndPoint * k1;

#if 1
    half yMax = csZBufferSize.y - 0.5;
    half yMin = 0.5;
    half xMax = csZBufferSize.x - 0.5;
    half xMin = 0.5;
    half alpha = 0;

    if (P1.y > yMax || P1.y < yMin)
    {
        half yClip = (P1.y > yMax) ? yMax : yMin;
        half yAlpha = (P1.y - yClip) / (P1.y - P0.y);
        alpha = yAlpha;
    }
    if (P1.x > xMax || P1.x < xMin)
    {
        half xClip = (P1.x > xMax) ? xMax : xMin;
        half xAlpha = (P1.x - xClip) / (P1.x - P0.x);
        alpha = max(alpha, xAlpha);
    }

    P1 = lerp(P1, P0, alpha);
    k1 = lerp(k1, k0, alpha);
    Q1 = lerp(Q1, Q0, alpha);
#endif

    P1 = (distanceSquared(P0, P1) < 0.0001) ? P0 + half2(0.01, 0.01) : P1;
    half2 delta = P1 - P0;
    bool permute = false;

    if (abs(delta.x) < abs(delta.y))
    {
        permute = true;
        delta = delta.yx;
        P1 = P1.yx;
        P0 = P0.yx;
    }

    half stepDirection = sign(delta.x);
    half invdx = stepDirection / delta.x;
    half2 dP = half2(stepDirection, invdx * delta.y);
    half3 dQ = (Q1 - Q0) * invdx;
    half dk = (k1 - k0) * invdx;
    
    dP *= stepSize;
    dQ *= stepSize;
    dk *= stepSize;
    P0 += dP * jitter;
    Q0 += dQ * jitter;
    k0 += dk * jitter;

    half3 Q = Q0;
    half k = k0;
    half prevZMaxEstimate = csOrigin.z;
    stepCount = 0;
    half rayZMax = prevZMaxEstimate, rayZMin = prevZMaxEstimate;
    half sceneZ = 100000;
    half end = P1.x * stepDirection;
    bool intersecting = intersectsDepthBuffer(rayZMin, rayZMax, sceneZ, layerThickness);
    half2 P = P0;
    int originalStepCount = 0;

    bool traceBehind_Old = true;
    rayIterations(forntDepth, traceBehind_Old, traceBehind, P, stepDirection, end, originalStepCount, maxSteps, intersecting, sceneZ, dP, Q, dQ, k, dk, rayZMin, rayZMax, prevZMaxEstimate, permute, hitPixel, invSize, layerThickness);

    stepCount = originalStepCount;
    Q.xy += dQ.xy * stepCount;
    csHitPoint = Q * (1 / k);
    return intersecting;
}

inline half3 ReconstructCSPosition(half4 _MainTex_TexelSize, half4 _ProjInfo, half2 S, half z)
{
    half linEyeZ = -LinearEyeDepth(z);
    return half3((((S.xy * _MainTex_TexelSize.zw)) * _ProjInfo.xy + _ProjInfo.zw) * linEyeZ, linEyeZ);
}

inline half3 GetPosition(sampler2D depth, half4 _MainTex_TexelSize, half4 _ProjInfo, half2 ssP)
{
    half3 P;
    P.z = SAMPLE_DEPTH_TEXTURE(depth, ssP.xy).r;
    P = ReconstructCSPosition(_MainTex_TexelSize, _ProjInfo, half2(ssP), P.z);
    return P;
}

/////////////////////////////////////Hierarchical_Z Trace/////////////////////////////////////
float2 cell(float2 ray, float2 cell_count) {
	return floor(ray.xy * cell_count);
}

float2 cell_count(float level, float2 ScreenSize) {
	return ScreenSize / (level == 0 ? 1 : exp2(level));
}

float3 intersectDepth_Plane(float3 rayOrigin, float3 rayDir, float marchSize)
{
	return rayOrigin + rayDir * marchSize;
}

float3 intersect_cell_boundary(float3 rayOrigin, float3 rayDir, float2 cellIndex, float2 cellCount, float2 crossStep, float2 crossOffset)
{
	float2 index = cellIndex + crossStep;  
        index /= cellCount;  
        index += crossOffset;

	float2 delta = index - rayOrigin.xy; 
        delta /= rayDir.xy;

	float marchSize = min(delta.x, delta.y);
	return intersectDepth_Plane(rayOrigin, rayDir, marchSize);
}

bool crossed_cell_boundary(float2 cell_id_one, float2 cell_id_two) {
	return (int)cell_id_one.x != (int)cell_id_two.x || (int)cell_id_one.y != (int)cell_id_two.y;
}

float minimum_depth_plane(float2 ray, float2 cell_count, float level, Texture2D SceneDepth) {
	return SceneDepth.Load( int3( (ray * cell_count), level ) );
}

float3 Hierarchical_Z_Trace(int HiZ_Max_Level, int HiZ_Start_Level, int HiZ_Stop_Level, int NumSteps, float2 screenSize, float3 rayOrigin, float3 rayDir, Texture2D SceneDepth)
{
	float level = HiZ_Start_Level;
	float2 crossStep = float2(rayDir.x >= 0.0 ? 1.0 : -1.0, rayDir.y >= 0.0 ? 1.0 : -1.0);
	float2 crossOffset = float2(crossStep.xy * screenSize); 
    crossStep.xy = saturate(crossStep.xy);

	float3 ray = rayOrigin.xyz;

	float3 marchSize = rayDir.xyz / rayDir.z;
	float3 curr_RayOrigin = intersectDepth_Plane(rayOrigin, marchSize, -rayOrigin.z);
	float2 rayCell = cell(ray.xy, screenSize);

	ray = intersect_cell_boundary(curr_RayOrigin, marchSize, rayCell.xy, screenSize, crossStep.xy, crossOffset.xy);
        
    int iterations = 0.0; float3 tmpRay = 0.0;
    
    [loop]
	while(level >=  HiZ_Stop_Level && iterations < NumSteps)
	{
        float2 cellCount = cell_count(level, screenSize);
		float minZ = minimum_depth_plane(ray.xy, cellCount, level, SceneDepth);
		float2 oldCellIdx = cell(ray.xy, cellCount);

		tmpRay = intersectDepth_Plane(curr_RayOrigin, marchSize, max(ray.z, minZ));
		float2 newCellIdx = cell(tmpRay.xy, cellCount);

        [branch]
		if(crossed_cell_boundary(oldCellIdx, newCellIdx))
		{
		    tmpRay = intersect_cell_boundary(curr_RayOrigin, rayDir, oldCellIdx, cellCount.xy, crossStep.xy, crossOffset.xy);
			level = min(HiZ_Max_Level, level + 2.0);
		}

		ray.xyz = tmpRay.xyz;
		level--;
		iterations++;
	}
	return ray;
}

float GetMarchSize(float2 start,float2 end,float2 SamplerPos)
{
    float2 dir = abs(end - start);
    return length( float2( min(dir.x, SamplerPos.x), min(dir.y, SamplerPos.y) ) );
}

float4 Hierarchical_Z_Trace(int HiZ_Max_Level, int HiZ_Start_Level, int HiZ_Stop_Level, int NumSteps, float thickness, float2 RayCastSize, float3 rayStart, float3 rayDir, Texture2D SceneDepth, SamplerState SceneDepth_Sampler)
{
    float SamplerSize = GetMarchSize(rayStart.xy, rayStart.xy + rayDir.xy, RayCastSize);
    float3 samplePos = rayStart + rayDir * (SamplerSize);
    int level = HiZ_Start_Level; float mask = 0.0;

    UNITY_LOOP
    for (int i = 0; i < NumSteps; i++)
    {
        float2 currSamplerPos = RayCastSize * exp2(level + 1.0);
        float newSamplerSize = GetMarchSize(samplePos.xy, samplePos.xy + rayDir.xy, currSamplerPos);
        float3 newSamplePos = samplePos + rayDir * newSamplerSize;
        float sampleMinDepth = Texture2DSampleLevel(SceneDepth, SceneDepth_Sampler, newSamplePos.xy, level);

        UNITY_FLATTEN
        if (sampleMinDepth < newSamplePos.z) {
            level = min(HiZ_Max_Level, level + 2.0);
            samplePos = newSamplePos;
        } else {
            level--;
        }

        UNITY_BRANCH
        if (level < HiZ_Stop_Level) {
            float delta = (-LinearEyeDepth(sampleMinDepth)) - (-LinearEyeDepth(samplePos.z));
            mask = delta <= thickness && i > 0.0;
            return float4(samplePos, mask);
        }
    }
    return float4(samplePos, mask);
}