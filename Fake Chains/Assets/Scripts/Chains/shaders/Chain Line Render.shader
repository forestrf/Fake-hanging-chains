Shader "Chain/Line Render " {
	Properties {
		_Color ("Color", Color) = (1,1,1,1)
		_MainTex ("Albedo (RGB)", 2D) = "white" {}
	}
	SubShader {
		Tags {
			"Queue" = "Geometry"
			"RenderType" = "Opaque"
		}
		LOD 200

		// ---- forward rendering base pass:
		Pass {
			Name "FORWARD"
			Tags { "LightMode" = "ForwardBase" }

			CGPROGRAM
			#pragma vertex vert_surf
			#pragma fragment frag_surf
			#pragma multi_compile_fwdbase
			
			#define UNITY_PASS_FORWARDBASE
			
			#include "HLSLSupport.cginc"
			#include "UnityShaderVariables.cginc"
			#include "UnityGlobalIllumination.cginc"
			#include "AutoLight.cginc"

			fixed4 _Color;
			sampler2D _MainTex; float4 _MainTex_ST;


			// vertex-to-fragment interpolation data
			struct v2f_surf {
				float4 pos : SV_POSITION;
				fixed4 color : COLOR;
				float2 pack0 : TEXCOORD0; // _MainTex
				half3 worldNormal : TEXCOORD1;
				fixed3 lightDir : TEXCOORD2;
				fixed3 worldPos : TEXCOORD3;
				
				SHADOW_COORDS(4)
			};

			// vertex shader
			v2f_surf vert_surf (appdata_full v) {
				v2f_surf o;
				UNITY_INITIALIZE_OUTPUT(v2f_surf, o);
				o.pos = mul (UNITY_MATRIX_MVP, v.vertex);
				o.color = v.color;
				
				o.pack0.xy = TRANSFORM_TEX(v.texcoord, _MainTex);
				
				o.worldPos = mul(_Object2World, v.vertex).xyz;
				o.worldNormal = UnityObjectToWorldNormal(v.normal);

				#ifndef USING_DIRECTIONAL_LIGHT
					o.lightDir = normalize(UnityWorldSpaceLightDir(o.worldPos));
				#else
					o.lightDir = _WorldSpaceLightPos0.xyz;
				#endif

				TRANSFER_SHADOW(o); // pass shadow coordinates to pixel shader
				
				return o;
			}

			// fragment shader
			fixed4 frag_surf (v2f_surf IN) : SV_Target {
				fixed4 c = 0;
				
				float2 uv_MainTex = IN.pack0.xy;

				UNITY_LIGHT_ATTENUATION(atten, IN, IN.worldPos)

				
				float3 Albedo = tex2D (_MainTex, uv_MainTex).rgb;
				Albedo *= _Color * IN.color;


				//#if (SHADER_TARGET < 30) || defined(SHADER_API_MOBILE)
				#ifdef UNITY_PASS_FORWARDBASE
					c.rgb = (_LightColor0.rgb * atten + UNITY_LIGHTMODEL_AMBIENT) * Albedo;
				#else
					c.rgb = (_LightColor0.rgb * atten) * Albedo;
					//c.rgb = _LightColor0.rgb * diff * 2 * atten;
				#endif

				return c;
			}
			ENDCG
		}

		// ---- SHADOW CASTING
		Pass {
			Name "ShadowCaster"
			Tags { "LightMode" = "ShadowCaster" }
			
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma multi_compile_shadowcaster
			
			#include "UnityCG.cginc"

			struct v2fF {
				V2F_SHADOW_CASTER;
			};

			struct appdataF {
				float4 vertex : POSITION;
				float3 normal : NORMAL;
				float2 texcoord0 : TEXCOORD0;
			};

			v2fF vert(appdataF v) {
				v2fF o;
				TRANSFER_SHADOW_CASTER_NORMALOFFSET(o)
				return o;
			}

			float4 frag( v2fF i ) : SV_Target {
				SHADOW_CASTER_FRAGMENT(i)
			}

			ENDCG
		}
	}
}
