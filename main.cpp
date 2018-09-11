#include <bits/stdc++.h>

using namespace std;

template <class T>
struct ValueType
{
	typedef typename T::value_type Result;
};

template <class T>
struct ValueType<T*>
{
	typedef T Result;
};

template <class T, size_t N>
struct ValueType<T[N]>
{
	typedef T Result;
};

enum AlgType { AlgType_Max = -1, AlgType_Min = 1 };

template <AlgType Type, typename Array>
void hungarian_alg( Array const& matrix, size_t n, vector<size_t>* distribution )
{
	typedef typename ValueType<Array>::Result ValueType;
	ValueType const Infinity = std::numeric_limits<ValueType>::max();

	ValueType u[n + 1], v[n + 1];
	for (size_t k = 0; k <= n; ++k)
		u[k] = v[k] = ValueType();

	size_t p[n + 1], way[n + 1];
	memset( p, 0, sizeof( size_t ) * (n + 1) );
	memset( way, 0, sizeof( size_t ) * (n + 1) );

	for (size_t i = 1; i <= n; ++i)
	{
		ValueType minv[n + 1];
		for (size_t k = 0; k <= n; ++k)
			minv[k] = Infinity;

		bool used[n + 1];
		memset( used, false, sizeof( bool ) * (n + 1) );

		p[0]      = i;
		size_t j0 = 0;
		do
		{
			used[j0]           = true;
			size_t const i0    = p[j0];
			size_t const row   = (i0 - 1) * n - 1;
			ValueType    delta = Infinity;
			size_t       j1    = 0;
			for (size_t j = 1; j <= n; ++j)
			{
				if (!used[j])
				{
					ValueType const cur = (Type == AlgType_Min) ? matrix[row + j] - u[i0] - v[j] : -matrix[row + j] - u[i0] - v[j];
					if (cur < minv[j])
					{
						minv[j] = cur;
						way[j]  = j0;
					}
					if (minv[j] < delta)
					{
						delta = minv[j];
						j1    = j;
					}
				}
			}
			for (size_t j = 0; j <= n; ++j)
			{
				if (used[j])
				{
					u[p[j]] += delta;
					v[j]    -= delta;
				}
				else
				{
					minv[j] -= delta;
				}
			}
			j0 = j1;
		} while (p[j0] != 0);

		do
		{
			size_t const j1 = way[j0];
			p[j0] = p[j1];
			j0    = j1;
		} while (j0);
	}

	distribution->resize( n );
	for (size_t j = 1; j <= n; ++j)
	{
		size_t const pj = p[j];
		if (pj != 0)
			(*distribution)[pj - 1] = j;
	}
}


class PriorityInt
{
public:
	PriorityInt() : _value(), _priority() {}
	PriorityInt( int value, int priority ) : _value( value ), _priority( priority ) {}

	PriorityInt& operator +=( PriorityInt const& other )       { _value += other._value; _priority += other._priority; return *this; }
	PriorityInt& operator -=( PriorityInt const& other )       { _value -= other._value; _priority -= other._priority; return *this; }
	PriorityInt  operator + ( PriorityInt const& other ) const { return PriorityInt( _value + other._value, _priority + other._priority ); }
	PriorityInt  operator - ( PriorityInt const& other ) const { return PriorityInt( _value - other._value, _priority - other._priority ); }
	PriorityInt  operator - () const                           { return PriorityInt( -_value, -_priority ); }

	bool         operator < ( PriorityInt const& other ) const { return _priority < other._priority || (_priority == other._priority && _value < other._value); }

private:
	int _value;
	int _priority;
};

namespace std
{
	template <>
	class numeric_limits<PriorityInt>
	{
	public:
		static PriorityInt max() { return PriorityInt( numeric_limits<int>::max(), numeric_limits<int>::max() ); }
	};
}

int main()
{
	vector<size_t> distribution;

	vector<int> a1 = {
		10, 20, 30,
		30, 30, 30,
		30, 30, 20,
	};
	hungarian_alg<AlgType_Min>( a1, 3, &distribution );
	for (auto a : distribution)
		cout << a << " ";
	cout << endl;

	vector<int> a2 = {
		7, 3, 6, 9, 5,
		7, 5, 7, 5, 6,
		7, 6, 8, 8, 9,
		3, 1, 6, 5, 7,
		2, 4, 9, 9, 5,
	};
	hungarian_alg<AlgType_Max>( a2, 5, &distribution );
	for (auto a : distribution)
		cout << a << " ";
	cout << endl;

	PriorityInt const affinity_table[] = {
		{ 0, 0 }, { 1, 5 }, { 0, 0 },
		{ 4, 3 }, { 1, 3 }, { 0, 0 },
		{ 1, 4 }, { 4, 4 }, { 0, 0 },
	};
	hungarian_alg<AlgType_Max>( affinity_table, 3, &distribution );
	for (auto a : distribution)
		cout << a << " ";
	cout << endl;

	return 0;
}
