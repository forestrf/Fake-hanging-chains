#if UNITY_EDITOR
using System;
using System.Collections;
#endif
using UnityEngine;

[ExecuteInEditMode]
public class FalseChainCreator : MonoBehaviour {
	[Header("Start and end points")]
	public Transform start;
	public Transform end;

	[Header("Link")]
	public Transform linkPrefab; // Create as many links as linkJoints
	public Vector3 linkForward = new Vector3(0, -1, 0);
	public Vector3 linkRotationAxis = new Vector3(0, 1, 0);

	[Header("Number of links")]
	public int numberOfLinks = 15;
	public float linkSize = 1;

	[Header("Other settings")]
	[Tooltip("Bigger values will be slower")]
	public int quality = 30; // Number of samples to capture the curve deformation. More = slower
	[Tooltip("Always update the position and rotation of the chain links")]
	public bool updateWhenOffscreen = false;


	//public AnimationCurve testCX, testCY;


	static GeometryUtilityUser guu = new GeometryUtilityUser();

	Transform[] linksTransforms;
	Vector2[] linkPositions; // update every frame looking at QuadRope and the Y function
	QuadRope qr = new QuadRope();
	int linkJoints; // Border links included
	float scale;
	Quaternion linkRotationAxisRotation;
	float ropeLength;
	Vector2 lastTarget;
	Vector3 lastTrgRot;



	void Start() {
		linkRotationAxisRotation = Quaternion.AngleAxis(90, linkRotationAxis);
		ropeLength = 1;
		//testCX = QuadRope.curveX;
		//testCY = QuadRope.curveY;
	}
	
	void LateUpdate() {
		//QuadRope.curveX = testCX;
		//QuadRope.curveY = testCY;

		linkJoints = numberOfLinks + 1;
		if (linkJoints < 2) return;
		if (linkSize <= 0) linkSize = 0.001f;
		if (null == start || null == end) return;
		if (null == linkPrefab) return;

		scale = numberOfLinks * linkSize;

		Vector3 sPos = start.position;
		Vector3 fromTo = end.position - sPos;
		Vector3 fromToScaled = fromTo;
		Vector3 lScale = transform.lossyScale;
		fromToScaled.x /= lScale.x;
		fromToScaled.y /= lScale.y;
		fromToScaled.z /= lScale.z;
		
		transform.position = sPos;
		Vector2 xz = new Vector2(fromToScaled.x, fromToScaled.z);
		Vector2 target = new Vector2(xz.magnitude, fromToScaled.y) / scale;
		if (target != lastTarget) {
			if (!updateWhenOffscreen && !guu.IsVisibleByCamera(new Bounds(sPos + fromTo * 0.5f, lScale * scale))) return;

			lastTarget = target;
			qr.SetTarget(target);
			
			SetLocalLinkPositions();
			
			SetLinkPositions();
		}
		Vector3 trgRot = new Vector3(xz.x, 0, xz.y);
		if (trgRot != lastTrgRot) {
			lastTrgRot = trgRot;
			transform.rotation = Quaternion.FromToRotation(Vector3.right, trgRot);
		}
	}

	string AnimationCurveToString(AnimationCurve a) {
		string s = "";
		for (int i = 0; i < a.length; i++) {
			s += "new Keyframe(" + a[i].time +
				"f, " + a[i].value +
				"f, " + a[i].inTangent +
				"f, " + a[i].outTangent +
				"f), ";
		}
		return s;
	}

	void OnDrawGizmosSelected() {
		if (linkJoints < 2) return;
		if (null == start || null == end) return;
		if (null == linkPrefab) return;

		//Debug.Log("x: " + AnimationCurveToString(testCX));
		//Debug.Log("y: " + AnimationCurveToString(testCY));

		Vector3 p = transform.position;

		// Draw box
		Gizmos.color = Color.white;
		Gizmos.DrawLine(p + (Vector3) qr.v[0], p + (Vector3) qr.v[1]);
		Gizmos.DrawLine(p + (Vector3) qr.v[1], p + (Vector3) qr.v[2]);
		Gizmos.DrawLine(p + (Vector3) qr.v[2], p + (Vector3) qr.v[3]);
		Gizmos.DrawLine(p + (Vector3) qr.v[3], p + (Vector3) qr.v[0]);

		// Draw rope curve
		Gizmos.color = Color.red;
		Vector3 lastP = p;
		for (int i = 0; i < quality; i++) {
			float i01 = (float) i / (quality - 1);
			Vector3 newP = p + (Vector3) qr.At(i01);
			Gizmos.DrawLine(lastP, newP);
			lastP = newP;
		}

		// Draw rope link joints
		Gizmos.color = Color.green;
		for (int i = 0; i < linkPositions.Length; i++) {
			Gizmos.DrawWireSphere(p + (Vector3) linkPositions[i], 0.02f);
		}

		// Update in editor to previsualize the setting
		linkRotationAxisRotation = Quaternion.AngleAxis(90, linkRotationAxis);

		/*
		float step = 0.1f;
		for (float y = 1; y >= -1; y -= step) {
			for (float x = 0; x <= 1; x += step) {
				if (x > Mathf.Cos(y)) continue;
				Vector2 v = new Vector2(x, y);
				qr.SetTarget(v);
				float l = GetRopeUnitLength();
				float l2 = Mathf.Pow(l, xx);
				l2 = 1 - Mathf.Abs(1 - l2);
				Gizmos.color = new Color(l2, l2, l2, 1);
				Gizmos.DrawWireSphere(p + (Vector3) v, step * 0.5f);
			}
		}
		*/
	}
	//public float xx;

	float GetRopeUnitLength() {
		return qr.GetCurveRectified();
	}

	void SetLinkPositions() {
		Vector3 lPos3 = linkPositions[0];
		for (int i = 0; i < linksTransforms.Length; i++) {
			Transform link = linksTransforms[i];
			Vector3 nextlPos3 = linkPositions[i + 1];

			Quaternion r = Quaternion.FromToRotation(linkForward, lPos3 - nextlPos3);
			if ((i & 1) == 0) { // Rotate 90º odd/even
				r = r * linkRotationAxisRotation;
			}

			link.localPosition = lPos3 * scale;
			link.localRotation = r;

			lPos3 = nextlPos3;
		}
	}

	void SetLocalLinkPositions() {
		if (linkPositions == null || linkPositions.Length != linkJoints) {
			linkPositions = new Vector2[linkJoints];
		}
		if (linksTransforms == null || linksTransforms.Length != numberOfLinks) {
			while (transform.childCount > 0)
				DestroyImmediate(transform.GetChild(0).gameObject);
			linksTransforms = new Transform[numberOfLinks];
			for (int i = 0; i < linksTransforms.Length; i++) {
				linksTransforms[i] = Instantiate(linkPrefab);
				linksTransforms[i].parent = transform;
				linksTransforms[i].localScale = Vector3.one;
			}
		}

		int link = 0;
		Vector2 v = qr.At(0);
		linkPositions[link++] = v;
		
		ropeLength = GetRopeUnitLength();

		float xOffset = 1f / (quality - 1);
		float linkDistance = 1f / numberOfLinks * ropeLength; // ropeLength [0..1]
		float distanciaRecorrida = 0;
		float x = 0;

		// Positionate last = target (approx)
		while (link < linkPositions.Length - 1) {
			float xNew = x + xOffset;
			Vector2 vNew = qr.At(xNew);
			float distanciaSegmento = Vector2.Distance(vNew, v);
			if (distanciaRecorrida + distanciaSegmento < linkDistance) {
				x = xNew;
				v = vNew;
				distanciaRecorrida += distanciaSegmento;
			} else {
				float lerp = (linkDistance - distanciaRecorrida) / distanciaSegmento;

				x = x + (xNew - x) * lerp; // Mathf.Lerp(x, xNew, lerp)
				v.x = v.x + (vNew.x - v.x) * lerp; // Mathf.Lerp(v, vNew, lerp)
				v.y = v.y + (vNew.y - v.y) * lerp; // Mathf.Lerp(v, vNew, lerp)
				linkPositions[link++] = v;
				distanciaRecorrida = 0;
			}
		}

		linkPositions[link] = qr.v[3];
	}
	

	
	/////////////////////////
	// CLASS
	/////////////////////////

	class QuadRope {
		public float width, height;

		public Vector2[] v = new Vector2[4];

		public static AnimationCurve curveX = new AnimationCurve(
			new Keyframe[] { new Keyframe(0f, 1f, -0.13f, -0.13f), new Keyframe(0.7082354f, 0.8226197f, -3.353576f, -3.353576f), new Keyframe(1f, 0f, -1.96f, -1.96f) });
		public static AnimationCurve curveY = new AnimationCurve(
			new Keyframe[] { new Keyframe(0f, 0f, 1.57f, 1.57f), new Keyframe(0.2524837f, 0.4345678f, 0.4180385f, 0.4180385f), new Keyframe(1f, 0.5f, 0f, 0f) });

		/*
		 * Project cuadratic funcion on a deformed quad. Start and end of the chain on target and base

		perimeter = 2
		(n) = Array index

		(3) + target
			|\_
			|  \_0.5
		 0.5|    \_
			|      \_+ base (0)
		(2) |        |
			 \_      |
			   \_    |0.5
			 0.5 \_  |
				   \_| (1)
		*/

		// 5 indices
		public void SetTarget(Vector2 targetOffset) {
			if (v.Length != 4) {
				Debug.LogError("quadPoints tiene que tener length = 4");
				return;
			}

			float distanceTargetBase = targetOffset.magnitude;

			if (distanceTargetBase > 1f) {
				// The chain will be stretched
				distanceTargetBase = 1;
			}

			// Lerp (4, 2, Mathf.Sqrt(distanceTargetBase)) // Not good enough results
			float div = 2 + 2 * curveX.Evaluate(distanceTargetBase); // Manually tweaked to make the chain length approx 1
			float distanceYAxis = (2 - distanceTargetBase * 2) / div;
			distanceYAxis = curveY.Evaluate(distanceYAxis);

			v[0] = Vector2.zero;
			v[1] = Vector2.zero - new Vector2(0, distanceYAxis);

			v[3] = targetOffset;
			v[2] = targetOffset - new Vector2(0, distanceYAxis);
		}

		// http://math.stackexchange.com/questions/389174/how-do-you-find-the-distance-between-two-points-on-a-parabola
		public float GetCurveRectified() {
			/*
			function:
			float g = v[3].x;
			float h = v[1].y;
			float j = v[3].y;
			x = x * 2 - 1
			f(x) := ((x / g * 2 - 1) ^ 2) * -h + j * x / g + h
			∫(√(1 + f'(x)^2),x)

			js to do the derive-to-c# faster
			String.prototype.replaceAll = function(search, replacement) {
				var target = this;
				return target.split(search).join(replacement);
			}
			a.replaceAll("^","").replaceAll("LN", "Mathf.Log").replaceAll("·", "*").replaceAll("√", "Mathf.Sqrt").replaceAll("ABS", "Mathf.Abs");
			*/

			//Integrate(1) - Integrate(0);
			return Integrate(v[3].x) - Integrate(0);
		}

		float Integrate(float x) {
			float g = v[3].x;
			float h = v[1].y;
			float j = v[3].y;
			
			float g2 = g * g;
			float hx = h * x;

			float jAdd4h = 4 * h + j;
			float sqrt1 = Mathf.Sqrt(64 * h * hx * x - 16 * g * hx * jAdd4h + g2 * g2 + g2 * jAdd4h * jAdd4h);
			float hx8_gjAdd4h = hx * 8 - g * jAdd4h;
			float h16 = 16 * h;

			return (Mathf.Log(sqrt1 + hx8_gjAdd4h) * g2 + sqrt1 * hx8_gjAdd4h / g2) / h16;
		}

		internal Vector2 At(float i01) {
			float g = v[3].x;
			float h = v[1].y;
			float j = v[3].y;

			float x = i01 * g;

			float xg21 = x / g * 2 - 1;
			return new Vector2(x, xg21 * xg21 * -h + j * x / g + h);
		}
	}

	class GeometryUtilityUser {
		PlaneSWR[] planes = new PlaneSWR[6] { new PlaneSWR(), new PlaneSWR(), new PlaneSWR(), new PlaneSWR(), new PlaneSWR(), new PlaneSWR() };
		int numberOfCameras = 0;
		Camera[] cameras = new Camera[1];
		Matrix4x4[] camarasVP = new Matrix4x4[2];
		int lastFrameCount;

		public bool IsVisibleByCamera(Bounds bounds) {
			if (lastFrameCount != Time.frameCount) {
				lastFrameCount = Time.frameCount;
				numberOfCameras = Camera.allCamerasCount;
				if (numberOfCameras > cameras.Length) {
					cameras = new Camera[numberOfCameras];
				}
#if UNITY_EDITOR
				ArrayList scenes = UnityEditor.SceneView.sceneViews;
				numberOfCameras += scenes.Count;
#endif

				if (numberOfCameras > camarasVP.Length) {
					camarasVP = new Matrix4x4[numberOfCameras];
				}

				for (int i = 0, l = Camera.GetAllCameras(cameras); i < l ; i++) {
					Camera cam = cameras[i];
					camarasVP[i] = cam.projectionMatrix * cam.worldToCameraMatrix;
				}
#if UNITY_EDITOR
				for (int i = 0; i < scenes.Count; i++) {
					Camera cam = ((UnityEditor.SceneView) scenes[i]).camera;
					camarasVP[numberOfCameras - i - 1] = cam.projectionMatrix * cam.worldToCameraMatrix;
				}
#endif
			}

			for (int i = 0; i < numberOfCameras; i++) {
				CalculateFrustumPlanes(planes, camarasVP[i]); // cachear proyeccion además de las cámaras
				if (CheckFrustumVisibility(planes, bounds)) {
					return true;
				}
			}

			return false;
		}

		// http://www.cnblogs.com/bodong/p/4800018.html
		// worldToProjectionMatrix = camera.projectionMatrix * camera.worldToCameraMatrix
		// Planes: 0 = Left, 1 = Right, 2 = Down, 3 = Up, 4 = Near, 5 = Far
		static void CalculateFrustumPlanes(PlaneSWR[] OutPlanes, Matrix4x4 worldToProjectionMatrix) {
			float RootVector0 = worldToProjectionMatrix.m30;
			float RootVector1 = worldToProjectionMatrix.m31;
			float RootVector2 = worldToProjectionMatrix.m32;
			float RootVector3 = worldToProjectionMatrix.m33;
			
			float ComVector0 = worldToProjectionMatrix.m00;
			float ComVector1 = worldToProjectionMatrix.m01;
			float ComVector2 = worldToProjectionMatrix.m02;
			float ComVector3 = worldToProjectionMatrix.m03;

			OutPlanes[0].Set(ComVector0 + RootVector0, ComVector1 + RootVector1, ComVector2 + RootVector2, ComVector3 + RootVector3);
			OutPlanes[1].Set(-ComVector0 + RootVector0, -ComVector1 + RootVector1, -ComVector2 + RootVector2, -ComVector3 + RootVector3);

			ComVector0 = worldToProjectionMatrix.m10;
			ComVector1 = worldToProjectionMatrix.m11;
			ComVector2 = worldToProjectionMatrix.m12;
			ComVector3 = worldToProjectionMatrix.m13;

			OutPlanes[2].Set(ComVector0 + RootVector0, ComVector1 + RootVector1, ComVector2 + RootVector2, ComVector3 + RootVector3);
			OutPlanes[3].Set(-ComVector0 + RootVector0, -ComVector1 + RootVector1, -ComVector2 + RootVector2, -ComVector3 + RootVector3);

			ComVector0 = worldToProjectionMatrix.m20;
			ComVector1 = worldToProjectionMatrix.m21;
			ComVector2 = worldToProjectionMatrix.m22;
			ComVector3 = worldToProjectionMatrix.m23;

			OutPlanes[4].Set(ComVector0 + RootVector0, ComVector1 + RootVector1, ComVector2 + RootVector2, ComVector3 + RootVector3);
			OutPlanes[5].Set(-ComVector0 + RootVector0, -ComVector1 + RootVector1, -ComVector2 + RootVector2, -ComVector3 + RootVector3);
		}

		//http://stagpoint.com/forums/threads/improved-intersection-tests-for-unity.8/
		/// <summary>
		/// Check if a bound is visible from a frustum
		/// </summary>
		static bool CheckFrustumVisibility(PlaneSWR[] frustumPlanes, Bounds bounds) {
			Vector3 boundsCenter = bounds.center;
			Vector3 boundExtents = bounds.extents;
			for (int i = 0; i < frustumPlanes.Length; i++) {
				Vector3 normal = frustumPlanes[i].normal;
				Vector3 normalAbs = frustumPlanes[i].normalAbs;

				// Compute the projection interval radius of box b onto L(t) = b.c + t * p.n
				//float r = Vector3.Dot(plane.normalAbs, boundExtents);
				float r = normalAbs.x * boundExtents.x + normalAbs.y * boundExtents.y + normalAbs.z * boundExtents.z;

				// Compute distance of box center from plane
				//float distance = Vector3.Dot(plane.normal, boundsCenter) + plane.distance;
				float distance = normal.x * boundsCenter.x + normal.y * boundsCenter.y + normal.z * boundsCenter.z + frustumPlanes[i].distance;


				if (distance < -r) {
					// If the AABB lies behind *any* of the planes, there
					// is no point in continuing with the rest of the test
					return false;
				}
			}

			return true;
		}

		class PlaneSWR {
			public Vector3 normal;
			public Vector3 normalAbs;
			public float distance;
			public PlaneSWR() { }
			public void Set(float InA, float InB, float InC, float InDistance) {
				normal.x = InA;
				normal.y = InB;
				normal.z = InC;

				float InverseMagnitude = 1.0f / normal.magnitude;
				normal.x *= InverseMagnitude;
				normal.y *= InverseMagnitude;
				normal.z *= InverseMagnitude;

				distance = InDistance * InverseMagnitude;
				normalAbs.x = normal.x > 0 ? normal.x : -normal.x;
				normalAbs.y = normal.y > 0 ? normal.y : -normal.y;
				normalAbs.z = normal.z > 0 ? normal.z : -normal.z;
			}
		}
	}
}
