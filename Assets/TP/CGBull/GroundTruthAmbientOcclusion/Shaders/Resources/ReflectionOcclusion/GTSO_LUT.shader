Shader "GTAO/GroundTruthSpecularOcclusion_LookUpTable"
{

	CGINCLUDE
		#include "UnityCG.cginc"
		#include "UnityCustomRenderTexture.cginc"
		#include "SpecularOcclusionIntegrate.cginc"


		fixed4 frag_Integrate(v2f_customrendertexture i) : SV_TARGET
		{
#if 1
			float3 uvw = i.localTexcoord.xyz;
			
			float thetaRef = uvw.x * 3.14 * 0.5;
			float roughness = uvw.y;

			float split = floor(uvw.z * 32);
			float cellZ = (split + 0.5) / 32.0;
			float cellW = uvw.z * 32 - split;
			float alphaV = 3.14 * 0.5 * cellZ;
			float beta = 3.14 * cellW;


			float GTSO_LUT = IntegrateGTSO(alphaV, beta, roughness, thetaRef);
			return GTSO_LUT;
#else
			float2 uv = i.localTexcoord.xy;
			float  alphaV = uv.x * 3.14 * 0.5;
			float thetaRef = uv.y * 3.14 * 0.5;
			float GTSO_LUT = IntegrateGTSO(alphaV, thetaRef, 1, thetaRef);
			return GTSO_LUT;
#endif
		}
	ENDCG
	SubShader
	{
		ZTest Always
		Cull Off
		ZWrite Off

		Pass 
		{ 
			CGPROGRAM 
				#pragma vertex CustomRenderTextureVertexShader
				#pragma fragment frag_Integrate
			ENDCG 
		}


	}
}

