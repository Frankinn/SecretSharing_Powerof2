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
#define PRIME_MAX 251
#define TEST_TIME 1
#define GF_MAX_EXP 255
#define GF_LENGTH 8
#define GF_EXP 0x100  //2^k
#define GF_IRPOLY 0x11B  //4: 0x13, 8: 0x11B/0x11D, 16: 0x1100B, 32: 0x100400007, 64: 0x1000000000000001B
#define GF_PRIME_DEGREE 7
#define TEST_DEGREE 10
using namespace std;
//g^k = a (ex. x^3 + 1 == 101 == 5 == a)
//short table[GF_MAX_EXP + 1];  //g^k --> a
//short arcTable[GF_MAX_EXP + 1];  //a --> g^k
//short inverseTable[GF_MAX_EXP + 1];  //g^k * g^j = e = 1,  g^k --> g^j
//bool setTable = false;
bool setPrimeTable = false;
bool completePrimeTable = false;
//int needNum;
int table_prime[251];  //GF_prime inverse table
const bool mx[] = { 1, 1, 0, 1, 1, 0, 0, 0, 1 }; // m(x) = x^8 + x^4 + x^3 + x + 1
const bool px[] = { 1, 1, 0, 1, 1, 1, 1, 1 };  // p = 251

random_device rd;
default_random_engine gen = default_random_engine(rd());
uniform_int_distribution<int> dis(0, GF_MAX_EXP);

class GaloisField {
public:
	//friend class Polynomial;
	//bool useTable = USE_TABLE;
	//short val = 0;
	bool arry[20] = { 0 };
	int degree = 0;

	GaloisField() {

	}

	GaloisField(const int n) {
		if (!n)
			return;
		for (int i = 0; i < GF_LENGTH; ++i)
			if(arry[i] = (n >> i) & 1)
				degree = i;
	}

	GaloisField(const bool irr[]) {
		for (int i = 0; i <= GF_LENGTH; ++i)
			arry[i] = mx[i];
		degree = GF_LENGTH;
	}

	GaloisField(const int s[], const int &n) {
		for (int i = 0; i <= n; ++i)
			arry[i] = s[i];
		degree = n;
	}

	GaloisField(const GaloisField &g) {
		degree = g.degree;
		for (int i = 0; i <= degree; ++i)
			arry[i] = g.arry[i];
	}

	const GaloisField &operator=(const GaloisField &g) {
		degree = g.degree;
		for (int i = 0; i <= degree; ++i)
			arry[i] = g.arry[i];
		for (int i = degree + 1; i <= GF_LENGTH; ++i)
			arry[i] = 0;
		return *this;
	}

	const GaloisField &operator=(const int &n) {
		if (!n) {
			degree = 0;
			for (int i = 0; i <= GF_LENGTH; ++i)
				arry[i] = 0;
			return *this;
		}
		for (int i = 0; i < GF_LENGTH; ++i)
			if (arry[i] = (n >> i) & 1)
				degree = i;
		return *this;
	}

	bool operator==(const GaloisField &g) const {
		if (degree != g.degree)
			return false;
		for (int i = 0; i <= degree; ++i)
			if (arry[i] != g.arry[i])
				return false;
		return true;
	}
	/*
	bool operator<=(const GaloisField &g) const {
		return degree <= g.degree;
	}
	
	bool operator<(const GaloisField &g) const {
		return degree < g.degree;
	}
	*/
	GaloisField operator+(const GaloisField &g) const {
		GaloisField ans;
		for (int i = 0; i <= degree; ++i)
			ans.arry[i] = arry[i] ^ g.arry[i];
		for (int i = degree + 1; i <= g.degree; ++i)
			ans.arry[i] = g.arry[i];
		ans.degree = degree > g.degree ? degree : g.degree;
		if(degree == g.degree)
			ans.normal();
		return ans;
	}

	GaloisField operator-(const GaloisField &g) const {
		GaloisField ans;
		for (int i = 0; i <= degree; ++i)
			ans.arry[i] = arry[i] ^ g.arry[i];
		for (int i = degree + 1; i <= g.degree; ++i)
			ans.arry[i] = g.arry[i];
		ans.degree = degree > g.degree ? degree : g.degree;
		if (degree == g.degree)
			ans.normal();
		return ans;
	}

	GaloisField operator*(const GaloisField &g) const {
		if (isZero() || g.isZero())
			return GaloisField();
		GaloisField ans;
		ans.degree = degree + g.degree;
		for (int i = 0; i <= degree; ++i) {
			if (!arry[i]) continue;
			for (int j = 0; j <= g.degree; ++j)
				ans.arry[i + j] ^= arry[i] && g.arry[j];
		}
		ans.normal();
		return ans;
	}

	GaloisField operator/(const GaloisField &g) const {
		/*if (*this < g)
			return GaloisField();*/
		return (*this * g.inverse());

	}

	const GaloisField &operator-() const {
		return *this;
	}

	bool &operator[](const int &index) {
		return arry[index];
	}

	void normal() {
		if (degree < 8) {
			for (int i = degree; i >= 0; --i)
				if (arry[i]) {
					degree = i;
					return;
				}
			degree = 0;
			return;
		}
		for (int i = degree; i >= GF_LENGTH; --i)
			if (arry[i]) {  //mx[] = { 1, 1, 0, 1, 1, 0, 0, 0, 1 }
				/*for (int j = GF_LENGTH, k = i; j >= 0; --j, --k)
					arry[k] ^= mx[j];*/
				arry[i] = 0;
				arry[i - 4] ^= 1;
				arry[i - 5] ^= 1;
				arry[i - 7] ^= 1;
				arry[i - 8] ^= 1;
			}
		for (int i = 7; i >= 0; --i)
			if (arry[i]) {
				degree = i;
				return;
			}
		degree = 0;
	}

	bool isZero() const {
		if (degree)
			return false;
		return !arry[0];
	}

	bool isOne() const {
		if (degree)
			return false;
		return arry[0];
	}

	const GaloisField &operator+=(const GaloisField &g) {
		*this = *this + g;
		return *this;
	}

	const GaloisField &operator-=(const GaloisField &g) {
		*this = *this - g;
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
		//Extended Euclid
		GaloisField a1 = 1, a2 = 0, a3 = mx;  //a3 = p = m(x) = x^8 + x^4 + x^3 + x + 1
		GaloisField b1 = 0, b2 = 1, b3 = *this;  //b3 = x = *this
		return extendedEuclid(a1, a2, a3, b1, b2, b3);
	}

	const GaloisField extendedEuclid(GaloisField a1, GaloisField a2, GaloisField a3, GaloisField b1, GaloisField b2, GaloisField b3) const {
		if (b3.isZero())
			return GaloisField();  //no inverse
		else if (b3.isOne())
			return b2;
		GaloisField q = 1 << (a3.degree - b3.degree);
		if(a3.degree >= 8)
			a3.normal();
		return extendedEuclid(b1, b2, b3, a1 - q * b1, a2 - q * b2, a3 - q * b3);
	}

	void newDegree() {
		for (int i = degree; i >= 0; --i)
			if (arry[i]) {
				degree = i;
				return;
			}
		degree = 0;
	}

	int val() const {
		int val = arry[degree];
		for (int i = degree - 1; i >= 0; --i) {
			val <<= 1;
			val |= arry[i];
		}
		return val;
	}

};

class GF_prime {
public:
	//int val = 0;
	//bool useTable = USE_TABLE;
	bool arry[20] = { 0 };
	int degree = 0;

	GF_prime() {
	}

	GF_prime(const int &n) {
		if (!n)
			return;
		for (int i = 0; i <= GF_PRIME_DEGREE; ++i)
			if (arry[i] = (n >> i) & 1)
				degree = i;
	}

	GF_prime(const bool p[]) {
		for (int i = 0; i <= GF_PRIME_DEGREE; ++i)
			arry[i] = px[i];
		degree = GF_PRIME_DEGREE;
	}

	const GF_prime operator+(const GF_prime &g) const;
	const GF_prime operator-(const GF_prime &g) const;
	const GF_prime operator*(const GF_prime &g) const;

	const GF_prime operator/(const GF_prime &g) const {
		return *this * g.inverse();
	}

	const GF_prime operator-() const;

	bool &operator[](const int &index) {
		if (index > degree)
			return arry[19];
		return arry[index];
	}

	const bool &operator<(const GF_prime &g) const {
		if (degree != g.degree)
			return degree < g.degree;
		for (int i = degree; i >= 0; --i) {
			if (arry[i] && !g.arry[i])  //*this > g
				return false;
			if (!arry[i] && g.arry[i])  //*this < g
				return true;
		}
		return false;
	}

	const bool &operator<=(const GF_prime &g) const {
		if (degree != g.degree)
			return degree < g.degree;
		for (int i = degree; i >= 0; --i) {
			if (arry[i] && !g.arry[i])
				return false;
			if (!arry[i] && g.arry[i])
				return true;
		}
		return true;
	}

	const bool &operator==(const GF_prime &g) const {
		if (degree != g.degree)
			return false;
		for (int i = 0; i <= degree; ++i)
			if (arry[i] != g.arry[i])
				return false;
		return true;
	}

	const GF_prime &operator=(const GF_prime &g) {
		degree = g.degree;
		for (int i = 0; i <= degree; ++i)
			arry[i] = g.arry[i];
		for (int i = degree + 1; i <= 16; ++i)
			arry[i] = 0;
		return *this;
	}

	const GF_prime &operator=(const int &n) {
		if (!n) {
			degree = 0;
			for (int i = 0; i <= 16; ++i)
				arry[i] = 0;
			return *this;
		}
		for (int i = 0; i <= GF_PRIME_DEGREE; ++i)
			if (arry[i] = (n >> i) & 1)
				degree = i;
		for (int i = GF_PRIME_DEGREE + 1; i <= 16; ++i)
			arry[i] = 0;
		return *this;
	}

	const GF_prime &operator+=(const GF_prime &g) {
		*this = *this + g;
		return *this;
	}

	const GF_prime &operator-=(const GF_prime &g) {
		*this = *this - g;
		return *this;
	}

	const GF_prime &operator*=(const GF_prime &g) {
		*this = *this * g;
		return *this;
	}

	const GF_prime &operator/=(const GF_prime &g) {
		*this = *this / g;
		return *this;
	}

	const GF_prime inverse() const {
		//Extended Euclid
		GF_prime a1 = 1, a2 = 0, a3 = px;
		GF_prime b1 = 0, b2 = 1, b3 = *this;
		return extendedEuclid(a1, a2, a3, b1, b2, b3);
	}

	const GF_prime extendedEuclid(GF_prime a1, GF_prime a2, GF_prime a3, GF_prime b1, GF_prime b2, GF_prime b3) const {
		if (b3.isZero())  //no inverse
			return GF_prime(0);
		else if (b3.isOne())
			return b2;
		GF_prime q = 1 << (a3.degree - b3.degree);
		if (a3 < b3 * q) {
			q = 1 << (a3.degree - b3.degree - 1);
		}
		return extendedEuclid(b1, b2, b3, a1 - q * b1, a2 - q * b2, a3 - q * b3);
	}

	void normal();

	const GF_prime operator<<(const int &n) const {  //shift left
		GF_prime ans;
		ans.degree = degree + n;
		for (int i = 0, j = n; i <= degree; ++i, ++j)
			ans.arry[j] = arry[i];
		return ans;
	}

	const GF_prime operator>>(const int &n) const {  //shift right
		GF_prime ans;
		ans.degree = degree - n;
		for (int i = n, j = 0; i <= degree; ++i, ++j)
			ans.arry[j] = arry[i];
		return ans;
	}

	void newDegree() {
		for (int i = degree; i >= 0; --i)
			if (arry[i]) {
				degree = i;
				return;
			}
		degree = 0;
	}

	bool isZero() const {
		if (degree)
			return false;
		return !arry[0];
	}

	bool isOne() const {
		if (degree)
			return false;
		return arry[0];
	}

	int val() const {
		int val = arry[degree];
		for (int i = degree - 1; i >= 0; --i) {
			val <<= 1;
			val |= arry[i];
		}
		return val;
	}
};

const GF_prime gfPrime = px;

const GF_prime GF_prime::operator+(const GF_prime &g) const {
	if (g.degree > degree)
		return g + *this;
	GF_prime ans;
	//ans.degree = degree > g.degree ? degree : g.degree;
	ans.degree = degree;
	int carry = 0;
	for (int i = 0; i <= g.degree; ++i) {
		carry += arry[i] + g.arry[i];
		ans[i] = carry % 2;
		carry /= 2;
	}
	for (int i = g.degree + 1; i <= degree; ++i) {
		carry += arry[i];
		ans[i] = carry % 2;
		carry /= 2;
	}
	if (carry) {
		++ans.degree;
		ans.arry[ans.degree] = 1;
	}
	/*if (gfPrime <= ans)
		ans -= gfPrime;*/
	ans.normal();
	return ans;
}

const GF_prime GF_prime::operator-(const GF_prime &g) const {
	if (isZero())
		return (-g);
	if (*this < g)
		return *this + (-g);
	GF_prime ans;
	//ans.degree = degree > g.degree ? degree : g.degree;
	ans.degree = degree;
	bool borrow = 0;
	for (int i = 0; i <= ans.degree; ++i) {
		if (arry[i] < g.arry[i] + borrow) {
			ans.arry[i] = arry[i] + 2 - g.arry[i] - borrow;
			borrow = 1;
		}
		else {
			ans.arry[i] = arry[i] - g.arry[i] - borrow;
			borrow = 0;
		}
	}
	ans.newDegree();
	return ans;
}

const GF_prime GF_prime::operator*(const GF_prime &g) const {
	if (isZero() || g.isZero())
		return GF_prime();
	if (isOne())
		return g;
	if (g.isOne())
		return *this;
	GF_prime ans;
	ans.degree = degree + g.degree;
	int carry = 0;
	for (int i = 0; i <= degree; ++i) {
		if (!arry[i]) continue;
		carry = 0;
		for (int j = 0; j <= g.degree; ++j) {
			carry += ans.arry[i + j] + g.arry[j];
			ans.arry[i + j] = carry % 2;
			carry /= 2;
		}
		if (carry >= 1) {
			ans.arry[i + g.degree + 1] = 1;
			if (i == degree)
				++ans.degree;
		}
	}
	ans.normal();
	return ans;
}

const GF_prime GF_prime::operator-() const {
	return gfPrime - *this;
}

void GF_prime::normal() {
	while (gfPrime < *this) {
		GF_prime p = gfPrime << (degree - GF_PRIME_DEGREE);
		if (*this < p)
			p = gfPrime << (degree - GF_PRIME_DEGREE - 1);
		*this -= p;
		newDegree();
	}
}

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

//int shareNum();

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
			/*if (poly[i] == 0) {
				--i;
				continue;
			}*/
		}
#if !TEST_MODE
		for (int i = degree; i >= 1; --i)
			cout << poly[i].val() << "x^" << i << " + ";  //print poly
		cout << poly[0].val() << endl;
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
		for (int i = degree; i >= 1; --i)
			cout << poly[i].val() << "x^" << i << " + ";  //print poly
		cout << poly[0].val() << endl;
#endif
	}

	void generateShare() {
		//generatePoly(secret);
		//Share *shares = new Share[degree + 1];  //degree + 1 shares
		int shareNum = peopleNum > degree + 1 ? peopleNum : degree + 1;
		shares = new Share<T>[shareNum];  //degree + 1 shares
		for (int i = 0; i < shareNum; ++i) {
			shares[i].x = dis(gen);
			if (checkShareRepeat<T>(shares, i) || shares[i].x.val() == 0) {
				--i;
				continue;
			}
			shares[i].y = poly[0];
			for (int j = 1; j <= degree; ++j) {
				shares[i].y = shares[i].y + poly[j] * gfPow<T>(shares[i].x, j);
			}
			if (shares[i].y.val() == 0) {
				--i;
				continue;
			}
		}
#if !TEST_MODE
		for (int i = 0; i < shareNum; ++i)
			cout << "(" << shares[i].x.val() << ' ' << shares[i].y.val() << ")  ";  //print shares
		cout << "\nfor " << peopleNum << " peoples, need " << degree + 1 << " shares to get " << secretNum << " secrets." << endl;
#endif
		//delete[] shares;
	}

	void printSecrets() {
		cout << "\nGet secrets: ";
		for (int i = secretNum - 1; i >= 0; --i)
			cout << poly[i].val() << ' ';
		cout << endl;
	}

	void printPoly() {
		//cout << "\n";
		int c = 97;
		for (int i = degree; i >= 1; --i, ++c)
			cout << poly[i].val() << ((char)c) << " + ";
		cout << poly[0].val() << ((char)c) << " = " << y.val() << endl;
	}
};

template<typename T>
int getSecret(const int &shareNum);  //Lagrange Interpolation, x = 0, get constant
template<typename T>
Polynomial< T > getAllSecret(const int &shareNum, const int &secretNum);  //Gaussian Elimination, get all coefficient

int main() {
	if (typeid(TYPE) == typeid(GaloisField))
		dis = uniform_int_distribution<int>(0, GF_MAX_EXP);
	else
		dis = uniform_int_distribution<int>(0, PRIME_MAX - 1);

	clock_t clock1, clock2;
	vector<int> clocks;
	Polynomial<TYPE> poly;
	double average = 0;
	///*
	int peopleNum = TEST_DEGREE + 1;
	int needNum = TEST_DEGREE + 1;
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
	cout << "type: " << TYPE_STR << "\nuse table: " << ((typeid(TYPE) == typeid(GaloisField)) ? (USE_TABLE ? "true" : "false") : (USE_TABLE ? "true" : "false")) << "\ndegree: " << needNum - 1 << "\ntest time: " << TEST_TIME << "\naverage time: " << average << "ms\n\n";
	clocks.clear();
	//*/

	
	//int secrets[] = { 20, 50, 11, 58, 9, 80, 105 };
	int secrets[] = { 9, 80, 105 };
	int secretNum = sizeof(secrets) / sizeof(secrets[0]);
	int peopleNum2 = TEST_DEGREE + 1;
	int needNum2 = TEST_DEGREE + 1;
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
	cout << "\ntype: " << TYPE_STR << "\nuse table: " << ((typeid(TYPE) == typeid(GaloisField)) ? (USE_TABLE ? "true" : "false") : (USE_TABLE ? "true" : "false")) << "\ndegree: " << needNum2 - 1 << "\ntest time: " << TEST_TIME << "\naverage time: " << average << "ms\n\n";
	clocks.clear();
	

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
	for (; i < exp; i *= 2)  //O(logn)
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
	return ans.val();
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
	for (int i = 0; i <= degree; ++i)  //test Gaussian
		poly[i].printPoly();
	*/
	for (int i = degree, d = degree; i > 0; --i, --d) {  //Gaussian Elimination, e.g. 4 poly, eliminate 3 others
		for (int j = i - 1; j >= 0; --j) {
			poly[j] *= poly[i].poly[d] / poly[j].poly[d];
			//poly[j] /= poly[j].poly[d];
			//poly[j] *= poly[i].poly[d];
			poly[j] -= poly[i];
		}

	}
	/*
	cout << endl;
	for (int i = 0; i <= degree; ++i)  //test Gaussian
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