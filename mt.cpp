#include"mt.h"

	 float Dot(const Vector3& v1, const Vector3& v2)
	{
		float result;
		result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

		return result;
	}

	 Vector3 Add(const Vector3& v1, const Vector3& v2)
	{
		Vector3 result{};
		result.x = v1.x + v2.x;
		result.y = v1.y + v2.y;
		result.z = v1.z + v2.z;
		return result;
	}

	 Vector3 Cross(const Vector3& v1, const Vector3& v2)
	{
		Vector3 result{};
		result.x = v1.y * v2.z - v1.z * v2.y;
		result.y = v1.z * v2.x - v1.x * v2.z;
		result.z = v1.x * v2.y - v1.y * v2.x;

		return result;
	}

	 Vector3 Subtract(const Vector3& v1, const Vector3& v2)
	{
		Vector3 result{};
		result.x = v1.x - v2.x;
		result.y = v1.y - v2.y;
		result.z = v1.z - v2.z;
		return result;
	}

	 Vector3 Multiply(const float& v1, const Vector3& v2)
	{
		Vector3 result{};
		result.x = v1 * v2.x;
		result.y = v1 * v2.y;
		result.z = v1 * v2.z;
		return result;
	}

	 Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		Matrix4x4 result{};

		for (int row = 0; row < 4; row++)
		{
			for (int column = 0; column < 4; column++)
			{
				result.m[row][column] = m1.m[row][0] * m2.m[0][column] + m1.m[row][1] * m2.m[1][column] + m1.m[row][2] * m2.m[2][column] + m1.m[row][3] * m2.m[3][column];

			}
		}
		return result;
	}

	 Matrix4x4 MakeTranslateMatrix(const Vector3 translate)
	{
		Matrix4x4 result{};

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result.m[i][j] = 0;
			}
		}

		result.m[0][0] = 1;
		result.m[1][1] = 1;
		result.m[2][2] = 1;
		result.m[3][3] = 1;


		result.m[3][0] = translate.x;
		result.m[3][1] = translate.y;
		result.m[3][2] = translate.z;

		return result;
	}

	 Matrix4x4 MakeScaleMatrix(const Vector3 scale)
	{
		Matrix4x4 result{};

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result.m[i][j] = 0;
			}
		}

		result.m[0][0] = scale.x;
		result.m[1][1] = scale.y;
		result.m[2][2] = scale.z;
		result.m[3][3] = 1;

		return result;
	}

	 Matrix4x4 MakeRotateXMatrix(float radian)
	{
		Matrix4x4 result{ };

		result.m[0][0] = 1;
		result.m[1][1] = std::cos(radian);
		result.m[1][2] = std::sin(radian);
		result.m[2][1] = -(std::sin(radian));
		result.m[2][2] = std::cos(radian);
		result.m[3][3] = 1;


		return result;
	}

	 Matrix4x4 MakeRotateYMatrix(float radian)
	{
		Matrix4x4 result{ };

		result.m[0][0] = std::cos(radian);
		result.m[0][2] = -(std::sin(radian));
		result.m[1][1] = 1;
		result.m[2][0] = std::sin(radian);
		result.m[2][2] = std::cos(radian);
		result.m[3][3] = 1;


		return result;
	}

	 Matrix4x4 MakeRotateZMatrix(float radian)
	{
		Matrix4x4 result{ };

		result.m[0][0] = std::cos(radian);
		result.m[1][0] = -(std::sin(radian));
		result.m[0][1] = std::sin(radian);
		result.m[1][1] = std::cos(radian);
		result.m[2][2] = 1;
		result.m[3][3] = 1;


		return result;
	}

	 Matrix4x4 MakeAffineMatrix(const Vector3 scale, const Vector3 rotate, const Vector3 translate)
	{
		Matrix4x4 result{};

		Matrix4x4 rotateXYZMatrix = Multiply(MakeRotateXMatrix(rotate.x), Multiply(MakeRotateYMatrix(rotate.y), MakeRotateZMatrix(rotate.z)));


		result.m[0][0] = rotateXYZMatrix.m[0][0] * scale.x;
		result.m[0][1] = rotateXYZMatrix.m[0][1] * scale.x;
		result.m[0][2] = rotateXYZMatrix.m[0][2] * scale.x;
		result.m[0][3] = 0;

		result.m[1][0] = rotateXYZMatrix.m[1][0] * scale.y;
		result.m[1][1] = rotateXYZMatrix.m[1][1] * scale.y;
		result.m[1][2] = rotateXYZMatrix.m[1][2] * scale.y;
		result.m[1][3] = 0;

		result.m[2][0] = rotateXYZMatrix.m[2][0] * scale.z;
		result.m[2][1] = rotateXYZMatrix.m[2][1] * scale.z;
		result.m[2][2] = rotateXYZMatrix.m[2][2] * scale.z;
		result.m[2][3] = 0;

		result.m[3][0] = translate.x;
		result.m[3][1] = translate.y;
		result.m[3][2] = translate.z;
		result.m[3][3] = 1;

		return result;
	}

	 Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRetio, float nearClip, float farClip)
	{
		Matrix4x4 result{};

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result.m[i][j] = 0;
			}
		}

		result.m[0][0] = (1 / aspectRetio) * (1 / tan(fovY / 2));
		result.m[1][1] = (1 / tan(fovY / 2));
		result.m[2][2] = farClip / (farClip - nearClip);
		result.m[2][3] = 1;
		result.m[3][2] = (-nearClip * farClip) / (farClip - nearClip);

		return result;
	}

	 Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip)
	{
		Matrix4x4 result{};

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result.m[i][j] = 0;
			}
		}

		result.m[0][0] = 2 / (right - left);
		result.m[1][1] = 2 / (top - bottom);
		result.m[2][2] = 1 / (farClip - nearClip);
		result.m[3][3] = 1;

		result.m[3][0] = (left + right) / (left - right);
		result.m[3][1] = (top + bottom) / (bottom - top);
		result.m[3][2] = (nearClip) / (nearClip - farClip);



		return result;
	}

	 Matrix4x4 MakeViewPortMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
	{
		Matrix4x4 result{};

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result.m[i][j] = 0;
			}
		}

		result.m[0][0] = width / 2;
		result.m[1][1] = -(height / 2);
		result.m[2][2] = maxDepth - minDepth;
		result.m[3][0] = left + (width / 2);
		result.m[3][1] = top + (height / 2);
		result.m[3][2] = minDepth;
		result.m[3][3] = 1;

		return result;
	}

	 Matrix4x4 Inverse(const Matrix4x4& m)
	{
		float a =
			m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3]
			+ m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1]
			+ m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]

			- m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1]
			- m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3]
			- m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]

			- m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3]
			- m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1]
			- m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]

			+ m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1]
			+ m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3]
			+ m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]

			+ m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3]
			+ m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1]
			+ m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]

			- m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1]
			- m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3]
			- m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]


			- m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0]
			- m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0]
			- m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]

			+ m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0]
			+ m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0]
			+ m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];


		Matrix4x4 result{};

		result.m[0][0] =
			(m.m[1][1] * m.m[2][2] * m.m[3][3]
				+ m.m[1][2] * m.m[2][3] * m.m[3][1]
				+ m.m[1][3] * m.m[2][1] * m.m[3][2]

				- m.m[1][3] * m.m[2][2] * m.m[3][1]
				- m.m[1][2] * m.m[2][1] * m.m[3][3]
				- m.m[1][1] * m.m[2][3] * m.m[3][2]
				) / a;

		result.m[0][1] =
			(-m.m[0][1] * m.m[2][2] * m.m[3][3]
				- m.m[0][2] * m.m[2][3] * m.m[3][1]
				- m.m[0][3] * m.m[2][1] * m.m[3][2]

				+ m.m[0][3] * m.m[2][2] * m.m[3][1]
				+ m.m[0][2] * m.m[2][1] * m.m[3][3]
				+ m.m[0][1] * m.m[2][3] * m.m[3][2]
				) / a;

		result.m[0][2] =
			(m.m[0][1] * m.m[1][2] * m.m[3][3]
				+ m.m[0][2] * m.m[1][3] * m.m[3][1]
				+ m.m[0][3] * m.m[1][1] * m.m[3][2]

				- m.m[0][3] * m.m[1][2] * m.m[3][1]
				- m.m[0][2] * m.m[1][1] * m.m[3][3]
				- m.m[0][1] * m.m[1][3] * m.m[3][2]
				) / a;

		result.m[0][3] =
			(-m.m[0][1] * m.m[1][2] * m.m[2][3]
				- m.m[0][2] * m.m[1][3] * m.m[2][1]
				- m.m[0][3] * m.m[1][1] * m.m[2][2]

				+ m.m[0][3] * m.m[1][2] * m.m[2][1]
				+ m.m[0][2] * m.m[1][1] * m.m[2][3]
				+ m.m[0][1] * m.m[1][3] * m.m[2][2]
				) / a;



		result.m[1][0] =
			(-m.m[1][0] * m.m[2][2] * m.m[3][3]
				- m.m[1][2] * m.m[2][3] * m.m[3][0]
				- m.m[1][3] * m.m[2][0] * m.m[3][2]

				+ m.m[1][3] * m.m[2][2] * m.m[3][0]
				+ m.m[1][2] * m.m[2][0] * m.m[3][3]
				+ m.m[1][0] * m.m[2][3] * m.m[3][2]
				) / a;

		result.m[1][1] =
			(m.m[0][0] * m.m[2][2] * m.m[3][3]
				+ m.m[0][2] * m.m[2][3] * m.m[3][0]
				+ m.m[0][3] * m.m[2][0] * m.m[3][2]

				- m.m[0][3] * m.m[2][2] * m.m[3][0]
				- m.m[0][2] * m.m[2][0] * m.m[3][3]
				- m.m[0][0] * m.m[2][3] * m.m[3][2]
				) / a;

		result.m[1][2] =
			(-m.m[0][0] * m.m[1][2] * m.m[3][3]
				- m.m[0][2] * m.m[1][3] * m.m[3][0]
				- m.m[0][3] * m.m[1][0] * m.m[3][2]

				+ m.m[0][3] * m.m[1][2] * m.m[3][0]
				+ m.m[0][2] * m.m[1][0] * m.m[3][3]
				+ m.m[0][0] * m.m[1][3] * m.m[3][2]
				) / a;

		result.m[1][3] =
			(m.m[0][0] * m.m[1][2] * m.m[2][3]
				+ m.m[0][2] * m.m[1][3] * m.m[2][0]
				+ m.m[0][3] * m.m[1][0] * m.m[2][2]

				- m.m[0][3] * m.m[1][2] * m.m[2][0]
				- m.m[0][2] * m.m[1][0] * m.m[2][3]
				- m.m[0][0] * m.m[1][3] * m.m[2][2]
				) / a;



		result.m[2][0] =
			(m.m[1][0] * m.m[2][1] * m.m[3][3]
				+ m.m[1][1] * m.m[2][3] * m.m[3][0]
				+ m.m[1][3] * m.m[2][0] * m.m[3][1]

				- m.m[1][3] * m.m[2][1] * m.m[3][0]
				- m.m[1][1] * m.m[2][0] * m.m[3][3]
				- m.m[1][0] * m.m[2][3] * m.m[3][1]
				) / a;

		result.m[2][1] =
			(-m.m[0][0] * m.m[2][1] * m.m[3][3]
				- m.m[0][1] * m.m[2][3] * m.m[3][0]
				- m.m[0][3] * m.m[2][0] * m.m[3][1]

				+ m.m[0][3] * m.m[2][1] * m.m[3][0]
				+ m.m[0][1] * m.m[2][0] * m.m[3][3]
				+ m.m[0][0] * m.m[2][3] * m.m[3][1]
				) / a;

		result.m[2][2] =
			(m.m[0][0] * m.m[1][1] * m.m[3][3]
				+ m.m[0][1] * m.m[1][3] * m.m[3][0]
				+ m.m[0][3] * m.m[1][0] * m.m[3][1]

				- m.m[0][3] * m.m[1][1] * m.m[3][0]
				- m.m[0][1] * m.m[1][0] * m.m[3][3]
				- m.m[0][0] * m.m[1][3] * m.m[3][1]
				) / a;

		result.m[2][3] =
			(-m.m[0][0] * m.m[1][1] * m.m[2][3]
				- m.m[0][1] * m.m[1][3] * m.m[2][0]
				- m.m[0][3] * m.m[1][0] * m.m[2][1]

				+ m.m[0][3] * m.m[1][1] * m.m[2][0]
				+ m.m[0][1] * m.m[1][0] * m.m[2][3]
				+ m.m[0][0] * m.m[1][3] * m.m[2][1]
				) / a;



		result.m[3][0] =
			(-m.m[1][0] * m.m[2][1] * m.m[3][2]
				- m.m[1][1] * m.m[2][2] * m.m[3][0]
				- m.m[1][2] * m.m[2][0] * m.m[3][1]

				+ m.m[1][2] * m.m[2][1] * m.m[3][0]
				+ m.m[1][1] * m.m[2][0] * m.m[3][2]
				+ m.m[1][0] * m.m[2][2] * m.m[3][1]
				) / a;

		result.m[3][1] =
			(m.m[0][0] * m.m[2][1] * m.m[3][2]
				+ m.m[0][1] * m.m[2][2] * m.m[3][0]
				+ m.m[0][2] * m.m[2][0] * m.m[3][1]

				- m.m[0][2] * m.m[2][1] * m.m[3][0]
				- m.m[0][1] * m.m[2][0] * m.m[3][2]
				- m.m[0][0] * m.m[2][2] * m.m[3][1]
				) / a;


		result.m[3][2] =
			(-m.m[0][0] * m.m[1][1] * m.m[3][2]
				- m.m[0][1] * m.m[1][2] * m.m[3][0]
				- m.m[0][2] * m.m[1][0] * m.m[3][1]

				+ m.m[0][2] * m.m[1][1] * m.m[3][0]
				+ m.m[0][1] * m.m[1][0] * m.m[3][2]
				+ m.m[0][0] * m.m[1][2] * m.m[3][1]
				) / a;

		result.m[3][3] =
			(m.m[0][0] * m.m[1][1] * m.m[2][2]
				+ m.m[0][1] * m.m[1][2] * m.m[2][0]
				+ m.m[0][2] * m.m[1][0] * m.m[2][1]

				- m.m[0][2] * m.m[1][1] * m.m[2][0]
				- m.m[0][1] * m.m[1][0] * m.m[2][2]
				- m.m[0][0] * m.m[1][2] * m.m[2][1]
				) / a;

		return result;
	}
	 Matrix4x4 Transpose(const Matrix4x4& m)
	{
		Matrix4x4 result{};

		for (int row = 0; row < 4; row++)
		{
			for (int column = 0; column < 4; column++)
			{
				result.m[row][column] = m.m[column][row];
			}
		}

		return result;
	}

	 Matrix4x4 MakeIdentity4x4()
	{
		Matrix4x4 result{};

		for (int row = 0; row < 4; row++)
		{
			for (int column = 0; column < 4; column++)
			{
				result.m[row][column] = 0;
			}
		}

		result.m[0][0] = 1;
		result.m[1][1] = 1;
		result.m[2][2] = 1;
		result.m[3][3] = 1;

		return result;
	}

	 Vector3 TransformCoord(Vector3 vector, Matrix4x4 matrix)
	{
		Vector3 result{};
		result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
		result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
		result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];

		float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
		assert(w != 0.0f);
		result.x /= w;
		result.y /= w;
		result.z /= w;
		return result;
	}

	 Vector3 Project(const Vector3& v1, const Vector3& v2)
	{
		Vector3 result{};

		result = Normalize(v2);
		float a = Dot(v1, result);
		result.x *= a;
		result.y *= a;

		return result;
	}

	 Vector3 ClosestPoint(const Vector3& point, const Segment& segment)
	{
		float t = Dot(Subtract(point, segment.origin), segment.diff) / std::powf(Length(segment.diff), 2.0f);
		Vector3 result = Add(segment.origin, Multiply(t, segment.diff));

		t = Clamp(t, 1.0f, 0.0f);


		return result;
	}

	 float Length(const Vector3& v)
	{
		float result;
		result = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);

		return result;
	}

	 Vector3 Normalize(const Vector3& v)
	{
		Vector3 result;
		result.x = v.x / Length(v);
		result.y = v.y / Length(v);
		result.z = v.z / Length(v);



		return result;
	}

	 float Clamp(float num, float max, float min)
	{
		if (num > max)
		{
			return max;
		}
		else if (num < min)
		{
			return min;
		}
		else
		{
			return num;
		}
	}
	 bool IsCollision(const Sphere& s1, const Sphere& s2)
	{
		float distance = Length(Subtract(s2.center, s1.center));

		if (distance <= s1.radius + s2.radius)
		{
			return	true;
		}

		return false;

	}

	 bool IsCollision(const Sphere& s1, const Plane& plane)
	{
		float k = std::abs(Dot(plane.normal, s1.center) - plane.distance);

		if (s1.radius > k)
		{
			return true;
		}

		return false;
	}


	 bool IsCollision2(const Segment& line, const Plane& plane)
	 {
		 float dot = Dot(plane.normal, line.diff);

		 if (dot == 0.0f)
		 {
			 return false;
		 }

		 float t = (plane.distance - Dot(line.origin, plane.normal)) / dot;

		 if (0.0f < t && t < 1.0f)
		 {
			 return true;
		 }

		 return false;
	 }



	 bool IsCollision3(const Triangle& triangle, const Segment& segment)
	 {
		 Plane plane{};
		 plane.normal =Normalize(Cross(Subtract(triangle.vertices[1], triangle.vertices[0]), Subtract(triangle.vertices[2], triangle.vertices[1])));

		 plane.distance = Dot(triangle.vertices[0], plane.normal);

		 float dot = Dot(plane.normal, segment.diff);

		 if (dot == 0.0f)
		 {
			 return false;
		 }

		 float t = (plane.distance - Dot(segment.origin, plane.normal)) / dot;

		 if (0.0f < t && t < 1.0f)
		 {
			 Vector3 p =Add(segment.origin,Multiply(t, segment.diff));

			 Vector3 cross01 = Cross(Subtract(triangle.vertices[1], triangle.vertices[0]), Subtract(p, triangle.vertices[1]));

			 Vector3 cross12 =Cross(Subtract(triangle.vertices[2], triangle.vertices[1]), Subtract(p, triangle.vertices[2]));

			 Vector3 cross20 = Cross(Subtract(triangle.vertices[0], triangle.vertices[2]), Subtract(p, triangle.vertices[0]));


			 if (Dot(cross01, plane.normal) >= 0.0f && Dot(cross12, plane.normal) >= 0.0f && Dot(cross20, plane.normal) >= 0.0f)
			 {
				 return true;
			 }
		 }

		 return false;
	 }


	 bool IsCollision4(const Over& aabb1, const Over& aabb2)
	 {
		 if ((aabb1.min.x <= aabb2.max.x && aabb1.max.x >= aabb2.min.x) &&
			 (aabb1.min.y <= aabb2.max.y && aabb1.max.y >= aabb2.min.y) &&
			 (aabb1.min.z <= aabb2.max.z && aabb1.max.z >= aabb2.min.z))
		 {

			 return true;

		 }

		 return false;
	 }

	 bool IsCollision5(const Over& aabb, const Sphere& sphere)
	 {
		 Vector3 clossestPoint
		 {
			 std::clamp(sphere.center.x, aabb.min.x, aabb.max.x),
			 std::clamp(sphere.center.y, aabb.min.y, aabb.max.y),
			 std::clamp(sphere.center.z, aabb.min.z, aabb.max.z)
		 };

		 float distance = Length(Subtract(clossestPoint, sphere.center));
		 if (distance <= sphere.radius)
		 {
			 return true;
		 }

		 return false;
	 }



	 Vector3 Perpendicular(const Vector3& vector)
	{
		if (vector.x != 0.0f || vector.y != 0.0f)
		{
			return { -vector.y, vector.x, 0.0f };
		}

		return { 0.0f, -vector.z, vector.y };
	}

	 void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
	{
		Vector3 center = Multiply(plane.distance, plane.normal);
		Vector3 perpendiculars[4];

		perpendiculars[0] = Normalize(Perpendicular(plane.normal));
		perpendiculars[1] = { -perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z };
		perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
		perpendiculars[3] = { -perpendiculars[2].x, -perpendiculars[2].y, -perpendiculars[2].z };

		Vector3 points[4];

		for (uint32_t index = 0; index < 4; ++index)
		{
			Vector3 extend = Multiply(2.0f, perpendiculars[index]);
			Vector3 point = Add(center, extend);
			points[index] = TransformCoord(TransformCoord(point, viewProjectionMatrix), viewportMatrix);
		}

		Novice::DrawLine(int(points[0].x), int(points[0].y), int(points[3].x), int(points[3].y), color);
		Novice::DrawLine(int(points[3].x), int(points[3].y), int(points[1].x), int(points[1].y), color);
		Novice::DrawLine(int(points[1].x), int(points[1].y), int(points[2].x), int(points[2].y), color);
		Novice::DrawLine(int(points[2].x), int(points[2].y), int(points[0].x), int(points[0].y), color);

	}

	 void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix)
	{
		const float kGridHalfWidth = 2.0f;
		const uint32_t kSubdivision = 10;
		const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);

		Vector3 localBorderVer[2]{};
		Vector3 localStripeVer[2]{};

		Vector3 screenBorderVer[2]{};
		Vector3 screenStripeVer[2]{};


		for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex)
		{


			localBorderVer[0] = { -kGridHalfWidth, 0.0f, kGridEvery * (float(xIndex) - 5) };
			localBorderVer[1] = { kGridHalfWidth, 0.0f, kGridEvery * (float(xIndex) - 5) };

			localStripeVer[0] = { kGridEvery * (int(xIndex) - 5) , 0.0f, -kGridHalfWidth };
			localStripeVer[1] = { kGridEvery * (int(xIndex) - 5) , 0.0f, kGridHalfWidth };


			Vector3 ndcBorderStart = TransformCoord(localBorderVer[0], viewProjectionMatrix);
			Vector3 ndcBorderEnd = TransformCoord(localBorderVer[1], viewProjectionMatrix);

			Vector3 ndcStripeStart = TransformCoord(localStripeVer[0], viewProjectionMatrix);
			Vector3 ndcStripeEnd = TransformCoord(localStripeVer[1], viewProjectionMatrix);


			screenBorderVer[0] = TransformCoord(ndcBorderStart, viewportMatrix);
			screenBorderVer[1] = TransformCoord(ndcBorderEnd, viewportMatrix);

			screenStripeVer[0] = TransformCoord(ndcStripeStart, viewportMatrix);
			screenStripeVer[1] = TransformCoord(ndcStripeEnd, viewportMatrix);

			Novice::DrawLine(int(screenBorderVer[0].x), int(screenBorderVer[0].y), int(screenBorderVer[1].x), int(screenBorderVer[1].y), 0xAAAAAAFF);

			Novice::DrawLine(int(screenStripeVer[0].x), int(screenStripeVer[0].y), int(screenStripeVer[1].x), int(screenStripeVer[1].y), 0xAAAAAAFF);

			if (localBorderVer[0].z == 0)
			{
				Novice::DrawLine(int(screenStripeVer[0].x), int(screenStripeVer[0].y), int(screenStripeVer[1].x), int(screenStripeVer[1].y), 0x000000FF);

				Novice::DrawLine(int(screenBorderVer[0].x), int(screenBorderVer[0].y), int(screenBorderVer[1].x), int(screenBorderVer[1].y), 0x000000FF);
			}

		}

	}

	 void DrawShere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
	{
		const uint32_t kSubdivision = 10;
		const float kLonEvery = (2 * 3.14f) / kSubdivision;
		const float kLatEvery = 3.14f / kSubdivision;

		for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex)
		{
			float lat = -3.14f / 2.0f + kLatEvery * latIndex;

			for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex)
			{
				float lon = lonIndex * kLonEvery;

				Vector3 a, b, c;
				a.x = sphere.center.x + sphere.radius * std::cosf(lat) * std::cosf(lon);
				a.y = sphere.center.y + sphere.radius * std::sinf(lat);
				a.z = sphere.center.z + sphere.radius * std::cosf(lat) * std::sinf(lon);

				b.x = sphere.center.x + sphere.radius * std::cosf(lat + (3.14f / kSubdivision)) * std::cosf(lon);
				b.y = sphere.center.y + sphere.radius * std::sinf(lat + (3.14f / kSubdivision));
				b.z = sphere.center.z + sphere.radius * std::cosf(lat + (3.14f / kSubdivision)) * std::sinf(lon);

				c.x = sphere.center.x + sphere.radius * std::cosf(lat) * std::cosf(lon + ((3.14f * 2) / kSubdivision));
				c.y = sphere.center.y + sphere.radius * std::sinf(lat);
				c.z = sphere.center.z + sphere.radius * std::cosf(lat) * std::sinf(lon + ((3.14f * 2) / kSubdivision));

				Vector3 ndcA = TransformCoord(a, viewProjectionMatrix);
				Vector3 ndcB = TransformCoord(b, viewProjectionMatrix);
				Vector3 ndcC = TransformCoord(c, viewProjectionMatrix);

				Vector3 screenA = TransformCoord(ndcA, viewportMatrix);
				Vector3 screenB = TransformCoord(ndcB, viewportMatrix);
				Vector3 screenC = TransformCoord(ndcC, viewportMatrix);

				Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
				Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenC.x), int(screenC.y), color);

			}
		}
	}

 void DrawLine(const Segment& seg, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
	 {
		 Vector3 start = TransformCoord(seg.origin, viewProjectionMatrix);
		 Vector3 screenStart = TransformCoord(start, viewportMatrix);
		 Vector3 end = TransformCoord(Add(seg.origin, seg.diff), viewProjectionMatrix);
		 Vector3 screenEnd = TransformCoord(end, viewportMatrix);
		 Novice::DrawLine(int(screenStart.x), int(screenStart.y), int(screenEnd.x), int(screenEnd.y), color);
	 }


 void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
 {
	 Vector3 vers[3]{};
	 Vector3 screenVers[3]{};

	 for (int i = 0; i < 3; i++)
	 {
		 vers[i] = TransformCoord(triangle.vertices[i], viewProjectionMatrix);

		 screenVers[i] = TransformCoord(vers[i], viewportMatrix);
	 }

	 Novice::DrawTriangle(int(screenVers[0].x), int(screenVers[0].y), int(screenVers[1].x), int(screenVers[1].y), int(screenVers[2].x), int(screenVers[2].y), color, kFillModeWireFrame);
 }

 void DrawOver(const Over& aabb, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
 {
	 Vector3 vers[8]{};
	 vers[0] = { aabb.min.x, aabb.min.y, aabb.min.z };
	 vers[1] = { aabb.min.x, aabb.min.y, aabb.max.z };
	 vers[2] = { aabb.min.x, aabb.max.y, aabb.min.z };
	 vers[3] = { aabb.max.x, aabb.min.y, aabb.min.z };
	 vers[4] = { aabb.max.x, aabb.max.y, aabb.min.z };
	 vers[5] = { aabb.min.x, aabb.max.y, aabb.max.z };
	 vers[6] = { aabb.max.x, aabb.min.y, aabb.max.z };
	 vers[7] = { aabb.max.x, aabb.max.y, aabb.max.z };

	 Vector3 screenVers[8]{};


	 for (int i = 0; i < 8; i++)
	 {
		 vers[i] =TransformCoord(vers[i], viewProjectionMatrix);
		 screenVers[i] =TransformCoord(vers[i], viewportMatrix);

	 }


	 Novice::DrawLine(int(screenVers[0].x), int(screenVers[0].y), int(screenVers[1].x), int(screenVers[1].y), color);
	 Novice::DrawLine(int(screenVers[0].x), int(screenVers[0].y), int(screenVers[2].x), int(screenVers[2].y), color);
	 Novice::DrawLine(int(screenVers[0].x), int(screenVers[0].y), int(screenVers[3].x), int(screenVers[3].y), color);

	 Novice::DrawLine(int(screenVers[1].x), int(screenVers[1].y), int(screenVers[5].x), int(screenVers[5].y), color);
	 Novice::DrawLine(int(screenVers[1].x), int(screenVers[1].y), int(screenVers[6].x), int(screenVers[6].y), color);

	 Novice::DrawLine(int(screenVers[2].x), int(screenVers[2].y), int(screenVers[4].x), int(screenVers[4].y), color);
	 Novice::DrawLine(int(screenVers[2].x), int(screenVers[2].y), int(screenVers[5].x), int(screenVers[5].y), color);

	 Novice::DrawLine(int(screenVers[3].x), int(screenVers[3].y), int(screenVers[4].x), int(screenVers[4].y), color);
	 Novice::DrawLine(int(screenVers[3].x), int(screenVers[3].y), int(screenVers[6].x), int(screenVers[6].y), color);

	 Novice::DrawLine(int(screenVers[4].x), int(screenVers[4].y), int(screenVers[7].x), int(screenVers[7].y), color);
	 Novice::DrawLine(int(screenVers[5].x), int(screenVers[5].y), int(screenVers[7].x), int(screenVers[7].y), color);
	 Novice::DrawLine(int(screenVers[6].x), int(screenVers[6].y), int(screenVers[7].x), int(screenVers[7].y), color);

 }
