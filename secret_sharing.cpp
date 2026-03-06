#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>
#include <typeinfo>
#define TEST_MODE false
#define TYPE GaloisField  //GaloisField GF_prime
#define TYPE_STR "GaloisField"
#define USE_TABLE false
#define TEST_TIME 1
#define PRIME_MAX 251
#define GF_LENGTH 8
#define GF_MAX_EXP 255
#define GF_EXP 0x100  //2^k
#define GF_IRPOLY 0x11B  //4: 0x13, 8: 0x11B/0x11D, 16: 0x1100B, 32: 0x100400007, 64: 0x1000000000000001B
using namespace std;
				 //g^k = a (ex. x^3 + 1 == 101 == 5 == a)
short table[GF_MAX_EXP+1];  //g^k --> a
short arcTable[GF_MAX_EXP+1];  //a --> g^k
short inverseTable[GF_MAX_EXP+1];  //g^k * g^j = e = 1,  g^k --> g^j
bool setTable = false;
bool setPrimeTable = false;
bool completePrimeTable = false;
//int needNum;
int table_prime[251];  //GF_prime inverse table

random_device rd;
default_random_engine gen = default_random_engine(rd());
uniform_int_distribution<int> dis(0, GF_MAX_EXP);

class GaloisField {
public:
	bool useTable = USE_TABLE;

	short val = 0;

	GaloisField() {
		if (useTable && !setTable)
			createTable();
		val = 0;
	}

	GaloisField(int n) {
		if (useTable && !setTable)
			createTable();
		val = n;
	}

	GaloisField(const GaloisField &g) {
		val = g.val;
	}

	const bool &operator==(const GaloisField &g) const {
		return (val == g.val);
	}

	const bool &operator!=(const GaloisField &g) const {
		return (val != g.val);  //return !(*this == g);
	}

	const GaloisField &operator=(const GaloisField &g) {
		val = g.val;
		return *this;
	}

	const GaloisField &operator=(const int &i) {
		val = i;
		return *this;
	}

	const GaloisField operator+(const GaloisField &g) const {
		return GaloisField(val ^ g.val);
	}

	const GaloisField operator-(const GaloisField &g) const {
		return GaloisField(val ^ g.val);
	}

	const GaloisField operator*(const GaloisField &g) const {
		if (!val || !g.val)
			return GaloisField();
		if(useTable)
			return GaloisField( table[ (arcTable[val] + arcTable[g.val]) % GF_MAX_EXP] );  //g^a * g^b == g^(a+b), val -> exp -> + -> val
		else {  //compute by bits
			GaloisField ans = 0, mul = g;
			for (int i = 0; i < GF_LENGTH; ++i) {
				if (mul.val & 1) {
					ans = ans + shiftLeft(*this, i);
				}
				mul.val = mul.val >> 1;
				if (mul.val == 0)
					return ans;
			}
			return ans;
		}
	}

	const GaloisField operator/(const GaloisField &g) const {
		return *this * g.inverse();  //division == multiply inverse, a/b = a * b^(-1)
	}

	const GaloisField &operator-() const {
		return *this;
	}

	const GaloisField &operator+=(const GaloisField &g) {
		val ^= g.val;
		return *this;
	}

	const GaloisField &operator-=(const GaloisField &g) {
		val ^= g.val;
		return *this;
	}

	const GaloisField &operator*=(const GaloisField &g) {
		*this = *this * g;
		return *this;
	}

	const GaloisField &operator/=(const GaloisField &g) {
		*this = *this / g;
		return *this;
	}

	const GaloisField inverse() const {
		if(useTable)
			return GaloisField(inverseTable[val]);
		else {  //Extended Euclid
			GaloisField a1 = 1, a2 = 0, a3 = GF_IRPOLY;  //a3 = p = m(x) = x^8 + x^4 + x^3 + x + 1
			GaloisField b1 = 0, b2 = 1, b3 = *this;  //b3 = x = *this
			return extendedEuclid(a1, a2, a3, b1, b2, b3);
		}
	}

	const GaloisField extendedEuclid(GaloisField a1, GaloisField a2, GaloisField a3, GaloisField b1, GaloisField b2, GaloisField b3) const {
		if (b3.val == 0)
			return GaloisField(-1);  //no inverse
		else if (b3.val == 1)
			return b2;
		GaloisField q = 1 << (a3.maxDegree() - b3.maxDegree());
		if (a3.val & GF_EXP)  //for first a3 = m(x) = x^8 + x^4 + x^3 + x + 1
			a3.val ^= GF_IRPOLY;  //only fist time a3 and q*b3 will bigger than max degree
		return extendedEuclid(b1, b2, b3, a1 - q*b1, a2 - q*b2, a3 - q*b3);
	}

	int maxDegree() const {
		int d = 0;
		GaloisField g = *this;
		for (int i = 0; i <= GF_LENGTH; ++i) {
			if (g.val & 1)
				d = i;
			g.val = g.val >> 1;
		}
		return d;
	}

	void createTable() {
		table[0] = 1;  //g^0
		for (int i = 1; i < GF_MAX_EXP; ++i) {  //g = x + 1
			table[i] = (table[i - 1] << 1) ^ table[i - 1];  //table[i-1] * (x + 1)
			if (table[i] & GF_EXP)  //max degree == 8
				table[i] ^= GF_IRPOLY;  //mod m(x) = x^8 + x^4 + x^3 + x + 1
		}

		for (int i = 0; i < GF_MAX_EXP; ++i)
			arcTable[table[i]] = i;

		for (int i = 1; i <= GF_MAX_EXP; ++i) {
			int k = (GF_MAX_EXP - arcTable[i]) % GF_MAX_EXP;
			inverseTable[i] = table[k];
		}
		setTable = true;
	}

	const GaloisField shiftLeft(const GaloisField &g, const int &num) const {
		GaloisField ans = g;
		for (int i = 0; i < num; ++i) {
			ans.val = ans.val << 1;
			if (ans.val & GF_EXP)
				ans.val ^= GF_IRPOLY;
		}
		return ans;
	}

	const GaloisField shiftRight(const GaloisField &g, const int &num) const {
		GaloisField ans = g;
		for (int i = 0; i < num; ++i) {
			ans.val = ans.val >> 1;
			/*if (ans.val == 0)
				return ans;*/
		}
		return ans;
	}

};
/*
class int_normal {
public:
	long long val = 0;
	int_normal(){}
	int_normal(const int &n) { val = n; }
	int_normal operator+(const int_normal &n) const { return int_normal(val + n.val); }
	int_normal operator-(const int_normal &n) const { return int_normal(val - n.val); }
	int_normal operator*(const int_normal &n) const { return int_normal(val * n.val); }
	int_normal operator/(const int_normal &n) const { return int_normal(val / n.val); }
	int_normal operator=(const int_normal &n) { return int_normal(val = n.val); }
	int_normal operator=(const int &n) { return int_normal(val + n); }
	int_normal operator+=(const int_normal &n) { val += n.val; return *this; }
	int_normal operator-=(const int_normal &n) { val -= n.val; return *this; }
	int_normal operator*=(const int_normal &n) { val *= n.val; return *this; }
	int_normal operator/=(const int_normal &n) { val /= n.val; return *this; }
	int_normal operator-() const { return int_normal(-val); }
	bool operator==(const int_normal &n) const { return val == n.val; }
};
*/
class GF_prime {
public:
	int val = 0;
	//int table[PRIME_MAX] = { 0 };  //inverse table
	bool useTable = USE_TABLE;
	GF_prime() {
		if (useTable && !setPrimeTable)
			createTable();
		val = 0;
	}

	GF_prime(const int &n) {
		if (useTable && !setPrimeTable)
			createTable();
		val = n;
	}

	const GF_prime operator+(const GF_prime &n) const {
		GF_prime ans = val + n.val;
		ans.val %= PRIME_MAX;
		return ans;
	}

	const GF_prime operator-(const GF_prime &n) const {
		GF_prime ans = val - n.val;
		while (ans.val < 0)
			ans.val += PRIME_MAX;
		ans.val %= PRIME_MAX;
		return ans;
	}

	const GF_prime operator*(const GF_prime &n) const {
		GF_prime ans = val * n.val;
		ans.val %= PRIME_MAX;
		return ans;
	}

	const GF_prime operator/(const GF_prime &n) const {
		GF_prime ans = *this * n.inverse();
		ans.val %= PRIME_MAX;
		return ans;
	}

	const GF_prime operator-() const {
		GF_prime ans = -val + PRIME_MAX;
		ans.val %= PRIME_MAX;
		return ans;
	}

	const bool &operator==(const GF_prime &n) const {
		return val == n.val;
	}

	const GF_prime &operator=(const GF_prime &n) {
		val = n.val;
		return *this;
	}

	const GF_prime &operator=(const int &n) {
		val = n;
		return *this;
	}

	const GF_prime &operator+=(const GF_prime &n) {
		*this = *this + n;
		return *this;
	}

	const GF_prime &operator-=(const GF_prime &n) {
		*this = *this - n;
		return *this;
	}

	const GF_prime &operator*=(const GF_prime &n) {
		*this = *this * n;
		return *this;
	}

	const GF_prime &operator/=(const GF_prime &n) {
		*this = *this / n;
		return *this;
	}

	const GF_prime inverse() const {
		if (useTable && completePrimeTable)
			return GF_prime(table_prime[val]);
		else {  //Extended Euclid
			GF_prime a1 = 1, a2 = 0, a3 = PRIME_MAX;
			GF_prime b1 = 0, b2 = 1, b3 = *this;
			return extendedEuclid(a1, a2, a3, b1, b2, b3);
		}
	}

	const GF_prime extendedEuclid(GF_prime a1, GF_prime a2, GF_prime a3, GF_prime b1, GF_prime b2, GF_prime b3) const {
		if (b3.val == 0)  //no inverse
			return GF_prime(0);
		else if (b3.val == 1)
			return b2;
		GF_prime q = a3.val / b3.val;
		return extendedEuclid(b1, b2, b3, a1 - q * b1, a2 - q * b2, a3 - q * b3);
	}

	void createTable() {
		setPrimeTable = true;
		for (int i = 0; i < PRIME_MAX; ++i)
			table_prime[i] = 0;
		for (int i = 1; i < PRIME_MAX; ++i) {
			if (!table_prime[i]) {
				int x = GF_prime(i).inverse().val;
				table_prime[i] = x;
				table_prime[x] = i;
			}
		}
		completePrimeTable = true;
	}
};

template<typename T>
struct Share {
	T x = 0;
	T y = 0;
};
Share<TYPE> *shares;  //temp

template<typename T>
T gfPow(const T &g, const int exp);
template<typename T>
bool checkShareRepeat(Share<T> *sh, int len);

template< typename T >
class Polynomial {
public:
	int degree = 0;
	int needNum = 0;
	int peopleNum = 0;
	int secretNum = 1;
	T *poly = nullptr;
	T y = 0;

	Polynomial() {

	}

	Polynomial(const int &d) {
		needNum = d;
		degree = d - 1;
		peopleNum = d;
		secretNum = 1;
		poly = new T[degree + 1];  //0 ~ degree
	}

	Polynomial(const int &need, const int &secrets, const int &people = 0) {
		needNum = need;
		secretNum = secrets;
		peopleNum = people;
		degree = need > secrets ? need - 1 : secrets - 1;
		poly = new T[degree + 1];
	}

	Polynomial(const Polynomial &p) {
		degree = p.degree;
		needNum = p.needNum;
		secretNum = p.secretNum;
		peopleNum = p.peopleNum;
		y = p.y;
		poly = new T[degree + 1];
		for (int i = 0; i <= degree; ++i)
			poly[i] = p.poly[i];
	}

	~Polynomial() {
		if (poly)
			delete[] poly;
		poly = nullptr;
	}

	T &operator[](const int &index) {
		return poly[index];
	}

	Polynomial operator+(const Polynomial &p) const {
		Polynomial ans = (*this);
		ans.y += p.y;
		for (int i = 0; i <= ans.degree; ++i)
			ans.poly[i] += p.poly[i];
		return ans;
	}

	Polynomial operator-(const Polynomial &p) const {
		Polynomial ans = (*this);
		ans.y -= p.y;
		for (int i = 0; i <= ans.degree; ++i)
			ans.poly[i] -= p.poly[i];
		return ans;
	}

	Polynomial operator*(const T &g) const {
		Polynomial ans = (*this);
		ans.y *= g;
		for (int i = 0; i <= ans.degree; ++i)
			ans.poly[i] *= g;
		return ans;
	}

	Polynomial operator/(const T &g) const {
		Polynomial ans = (*this);
		ans.y /= g;
		for (int i = 0; i <= ans.degree; ++i)
			ans.poly[i] /= g;
		return ans;
	}

	const Polynomial &operator+=(const Polynomial &p) {
		*this = *this + p;
		return (*this);
	}

	const Polynomial &operator-=(const Polynomial &p) {
		*this = *this - p;
		return (*this);
	}

	const Polynomial &operator*=(const T &g) {
		*this = *this * g;
		return (*this);
	}

	const Polynomial &operator/=(const T &g) {
		*this = *this / g;
		return (*this);
	}

	const bool operator==(const Polynomial &p) const {
		if (degree != p.degree || (needNum != p.needNum) || y != p.y)
			return false;
		for (int i = 0; i <= degree; ++i)
			if (poly[i] != p.poly[i])
				return false;
		return true;
	}

	const Polynomial &operator=(const Polynomial &p) {
		degree = p.degree;
		needNum = p.needNum;
		secretNum = p.secretNum;
		peopleNum = p.peopleNum;
		y = p.y;
		if (poly)
			delete[] poly;
			
		poly = new T[degree + 1];
		for (int i = 0; i <= degree; ++i)
			poly[i] = p.poly[i];
		return (*this);
	}

	void generatePoly(const T &secret) {  //one secret
		poly[0] = secret;
		for (int i = 1; i <= degree; ++i) {
			poly[i] = dis(gen);
		}
#if !TEST_MODE
		///*
		for (int i = degree; i >= 1; --i)
			cout << poly[i].val << "x^" << i << " + ";  //print poly
		cout << poly[0].val << endl;
		//*/
#endif
	}

	void generatePoly(const int secrets[], const int &secretNum) {  //many secrets
		//int needNum = degree + 1;
		if (needNum < secretNum) {
			degree = secretNum - 1;
			if (poly)
				delete[] poly;
			poly = new T[degree + 1];
		}

		for (int i = 0; i < secretNum; ++i)
			poly[i] = secrets[i];
		for (int i = secretNum; i <= degree; ++i) {
			poly[i] = dis(gen);
			/*if (poly[i] == 0) {
				--i;
				continue;
			}*/
		}
#if !TEST_MODE
		///*
		for (int i = degree; i >= 1; --i)
			cout << poly[i].val << "x^" << i << " + ";  //print poly
		cout << poly[0].val << endl;
		//*/
#endif
	}

	void generateShare() {
		//generatePoly(secret);
		//Share *shares = new Share[degree + 1];  //degree + 1 shares
		int shareNum = peopleNum > degree + 1 ? peopleNum : degree + 1;
		shares = new Share<T>[shareNum];  //degree + 1 shares
		for (int i = 0; i < shareNum; ++i) {
			shares[i].x = dis(gen);
			if (checkShareRepeat<T>(shares, i) || shares[i].x.val == 0) {
				--i;
				continue;
			}
			shares[i].y = poly[0];
			for (int j = 1; j <= degree; ++j) {
				shares[i].y = shares[i].y + poly[j] * gfPow<T>(shares[i].x, j);
			}
			if (shares[i].y.val == 0) {
				--i;
				continue;
			}
		}
#if !TEST_MODE
		///*
		for (int i = 0; i < shareNum; ++i)
			cout << "(" << shares[i].x.val << ' ' << shares[i].y.val << ")  ";  //print shares
		cout << "\nfor " << peopleNum << " peoples, need " << degree + 1 << " shares to get " << secretNum << " secrets." << endl;
		//*/
#endif
		//delete[] shares;
	}

	void printSecrets() {
		cout << "\nGet secrets: ";
		for (int i = secretNum - 1; i >= 0; --i)
			cout << poly[i].val << ' ';
		cout << endl;
	}

	void printPoly() {
		//cout << "\n";
		int c = 97;
		for (int i = degree; i >= 1; --i, ++c)
			cout << poly[i].val << ((char)c) << " + ";
		cout << poly[0].val << ((char)c) << " = " << y.val << endl;
	}
};

template<typename T>
int getSecret(const int &shareNum);  //Lagrange Interpolation, x = 0, get constant
template<typename T>
Polynomial< T > getAllSecret(const int &shareNum, const int &secretNum);  //Gaussian Elimination, get all coefficient

int main() {
	if(typeid(TYPE) == typeid(GaloisField))
		dis = uniform_int_distribution<int>(0, GF_MAX_EXP);
	else
		dis = uniform_int_distribution<int>(0, PRIME_MAX - 1);

	clock_t clock1, clock2;
	vector<int> clocks;
	Polynomial<TYPE> poly;
	double average = 0;
	///*
	int peopleNum = 11;
	int needNum = 11;
	int secret = 105;
	for (int i = 0; i < TEST_TIME; ++i) {
	clock1 = clock();
	poly = Polynomial<TYPE>(needNum, 1, peopleNum);
	poly.generatePoly(secret);
	poly.generateShare();
#if TEST_MODE
	getSecret<TYPE>(needNum);  //test temp
#else
	cout << "\nGet Secret: " << getSecret<TYPE>(needNum) << endl << endl;
#endif
	clock2 = clock();
	if (shares)
		delete[] shares;
	poly.~Polynomial();
	clocks.push_back(clock2 - clock1);
	//cout << "type: " << TYPE_STR << "\ndegree: " << poly.degree << "\ntime: " << clock2 - clock1 << "ms\n\n";
	}

	for (int i = 0; i < TEST_TIME; ++i)
		average += clocks[i];
	average /= TEST_TIME;
	//cout << clock2 << ' ' << clock1 << endl;
	cout << "type: " << TYPE_STR << "\nuse table: " << ((typeid(TYPE) == typeid(GaloisField)) ? (USE_TABLE ? "true" : "false") : (USE_TABLE ? "true" : "false")) << "\ndegree: " << needNum-1 << "\ntest time: " << TEST_TIME << "\naverage time: " << average << "ms\n\n";
	clocks.clear();
	//*/
	
	///*
	//int secrets[] = { 20, 50, 11, 58, 9, 80, 105 };
	int secrets[] = { 9, 80, 105 };
	int secretNum = sizeof(secrets) / sizeof(secrets[0]);
	int peopleNum2 = 11;
	int needNum2 = 11;
	Polynomial<TYPE> morePoly;
	for (int i = 0; i < TEST_TIME; ++i) {
	clock1 = clock();
	//vector<int> secrets(s);
	morePoly = Polynomial<TYPE>(needNum2, secretNum, peopleNum2);
	morePoly.generatePoly(secrets, secretNum);
	morePoly.generateShare();
	Polynomial<TYPE> secretsPoly = getAllSecret<TYPE>(needNum2, secretNum);
	clock2 = clock();
	if (shares)
		delete[] shares;
	morePoly.~Polynomial();
#if !TEST_MODE
	secretsPoly.printSecrets();
#endif
	clocks.push_back(clock2 - clock1);
	}

	average = 0;
	for (int i = 0; i < TEST_TIME; ++i)
		average += clocks[i];
	average /= TEST_TIME;

	//cout << clock2 << ' ' << clock1 << endl;
	cout << "\ntype: " << TYPE_STR << "\nuse table: " << ((typeid(TYPE) == typeid(GaloisField)) ? (USE_TABLE ? "true" : "false") : (USE_TABLE ? "true" : "false")) << "\ndegree: " << needNum2-1 << "\ntest time: " << TEST_TIME << "\naverage time: " << average << "ms\n\n";
	clocks.clear();
	//*/

	//system("pause");
}

template<typename T>
T gfPow(const T &g, const int exp) {
	if (!exp)
		return T(1);
	if (exp == 1)
		return g;
	T ans = g * g;
	int i = 4;
	for (; i < exp; i*=2)  //O(logn)
		ans = ans * ans;
	for (i /= 2; i < exp; ++i)
		ans *= g;
	return ans;
}

template<typename T>
bool checkShareRepeat(Share<T> *sh, int len) {
	for (int i = 0; i < len; ++i)
		if (sh[i].x == sh[len].x)
			return true;
	return false;
}

template<typename T>
int getSecret(const int &shareNum) {  //Lagrange Interpolation, x = 0, get constant coefficient
	T ans = 0;
	for (int i = 0; i < shareNum; ++i) {
		T li = 1, lu = 1, ld = 1;
		for (int j = 0; j < shareNum; ++j) {
			if (j == i) continue;
			lu = lu * (-shares[j].x);  // "-" for GF_prime
			ld = ld * (shares[i].x - shares[j].x);
		}
		li = shares[i].y * lu / ld;
		ans = ans + li;
	}
	return ans.val;
}

template<typename T>
Polynomial<T> getAllSecret(const int &needNum, const int &secretNum) {  //Gaussian Elimination, get all coefficient
	int degree = needNum > secretNum ? needNum - 1 : secretNum - 1;
	Polynomial<T> *poly = new Polynomial<T>[degree + 1];
	for (int i = 0; i <= degree; ++i) {  //input shares in to function
		poly[i] = Polynomial<T>(needNum, secretNum); //				 degree  3    2    1    0   y
		for (int j = 0; j <= degree; ++j)  //ax^3 + bx^2 + cx + d = y  -->  ?a + ?b + ?c + ?d = ?
			poly[i].poly[j] = gfPow<T>(shares[i].x, j);
		poly[i].y = shares[i].y;
	}
	/*
	cout << endl;
	for (int i = 0; i <= degree; ++i)  //test
		poly[i].printPoly();
	*/
	for (int i = degree, d = degree; i > 0; --i, --d) {  //Gaussian Elimination, e.g. 4 poly, eliminate 3 others
		for (int j = i - 1; j >= 0; --j) {
			poly[j] *= poly[i].poly[d] / poly[j].poly[d];
			poly[j] -= poly[i];
		}

	}
	/*
	cout << endl;
	for (int i = 0; i <= degree; ++i)  //test
		poly[i].printPoly();
	*/
	Polynomial<T> ans(needNum, secretNum);
	ans[0] = poly[0].y / poly[0].poly[0];  //last poly, e.g. d = y
	for (int i = 1; i <= degree; ++i) {  //d = y  -->  ?c + d = y  -->  ?b + ?c d = y  -->  ?a + ?b + ?c d = y
		ans[i] = poly[i].y;
		for (int j = 0; j < i; ++j)
			ans[i] -= ans[j] * poly[i].poly[j];
		ans[i] /= poly[i].poly[i];
	}
	return ans;
}