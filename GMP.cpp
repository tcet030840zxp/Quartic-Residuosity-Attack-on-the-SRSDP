// GMP.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <gmp.h>
#include <gmpxx.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include "gaussian_integers.h"

using namespace std;
/*int record[10000000] = {0};
int main()
{
	freopen("r.txt", "r", stdin);
	freopen("w.txt", "w", stdout);
	GI alpha{ 1, 2 };
	GI beta{ -3, 5 };
	cout << alpha + beta << endl;
	cout << alpha - beta << endl;
	cout << alpha * beta << endl;
	cout << beta / alpha << endl;
	cout << beta % alpha << endl;
	cout << GIgcd({ 347, 89 }, { 117, -547 }) << endl;
	vector<GI> v = GIexgcd({ 347, 89 }, { 117, -547 });
	for (GI a : v)
		cout << a << endl;
	cout << GImodlinearsolve({ 347, 89 }, { 117, -547 }, { 2003 }) << endl;
	GImodlinearsolve({ 347, 89 }, { 117, -547 }, { 1, 2 });
	printGIproduct(GIprimefactor({ 1284, -418764 }));
	cout << QuarticResSymbol({ 347, 89 }, { 5, 6 }) << endl;
	cout << isGIprime({ 5, 6 }) << endl;

	cout << "开始" << endl;
	//mpz p0 = 2, ps = 5, pt = 11, qs = 7, qt = 19; //GCD(g-i,N)=Norm(p);
	mpz p0 = 2, ps = 5, pt = 47, qs = 7, qt = 19;
	int d = 3;
	mpz p0d = pow(p0, d);
	mpz p = 2 * p0d * ps * pt + 1;
	mpz q = 2 * p0d * qs * qt + 1;
	mpz N = p * q;
	cout << "ps的二进制表示是" << decToBin(ps) << endl;
	cout << "qs的二进制表示是" << decToBin(qs) << endl;
	cout << "pt的二进制表示是" << decToBin(pt) << endl;
	cout << "qt的二进制表示是" << decToBin(qt) << endl;
	mpz phiN = Eulerphi(N);
	factor(N);
	factor(p);
	factor(q);
	cout << "N = ";
	printGIproduct(GIprimefactor({ N,0 }));
	cout << "Eulerphi(N) = " << phiN.get_str() << endl;
	mpz g = 1;
	for (; g < N; g++)
	{
		if (isorder(g, p0d, N))
			break;
	}
	cout << "g = " << g.get_str() << endl;
	g = modularExponentiation(g, p0d / 4, N);
	cout << "g = " << g.get_str() << endl;
	GI gamma = { g, -1 };
	cout << "GCD(" << N.get_str() << ", " << gamma << ") = ";
	beta = GIgcd(gamma, { N, 0 });
	printGIproduct(GIprimefactor(beta));
	int cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;// 1 -1 i -i
	for (int i = 1; i < N; i++)
	{
		GI alpha{ i, 0 };
		if (JacobiSymbol(i, p) == 1 && JacobiSymbol(i, q) == 1)
		{
			gamma = QuarticResSymbol(alpha, beta);
			if (gamma == GI(1, 0))cnt1++;
			else if (gamma == GI(-1, 0))
			{
				cnt2++;
				if (i > 200000)
				{
					for (int j = 1; j < N; j++)
					{
						if (i == modularExponentiation(j, 2, N))
						{
							cout << "sqrt(x) = " << j << endl;
							break;
						}
					}
					for (int j = 1;; j++)
					{
						if (modularExponentiation(i, j, N) == 1)
						{
							cout << "x的阶是" << j << endl;
							break;
						}
					}
					for (int j = 1; j < N; j++)
					{
						if (modularExponentiation(j, p0d * pt * qt, N) == i)
						{
							cout << "x有形式y^{p_0^dp_tq_t}" << j << endl;
							break;
						}
					}
					QuarticResSymbol(alpha, beta, 1);
					break;
				}
			}
			else if (gamma == GI(0, 1))cnt3++;
			else cnt4++;
		}
	}
	cout << "QR_N中四次剩余符号是1的数目: " << cnt1 << endl;
	cout << "QR_N中雅克比符号和四次剩余符号是-1的数目: " << cnt2 << endl;
	cout << "QR_N中雅克比符号和四次剩余符号是i的数目: " << cnt3 << endl;
	cout << "QR_N中雅克比符号和四次剩余符号是-i的数目: " << cnt4 << endl;
	int cnt5 = 0, cnt6 = 0;
	for (int i = 1; i < N; i++)
	{
		//if (i % 10000 == 0)cout << i << endl;
		if (JacobiSymbol(i, p) == 1 && JacobiSymbol(i, q) == 1)
		{
			mpz tmp = 1;
			for (int j = 1; j <= ps * qs; j++)
			{
				tmp = (tmp * i) % N;
				if (j == ps * qs && tmp == 1)cnt5++;
				if (j < ps * qs && tmp == 1)break;
			}
			record[modularExponentiation(i, p0d * pt * qt, N).get_ui()]++;
		}
	}
	cout << cnt5 << endl;
	for (int i = 1; i < N; i++)if (record[i])cnt6++;
	cout << cnt6 << endl;
	return 0;
}
const int sz = 7680;
int main()
{
	//生成大素数算法
	gmp_randstate_t grt;
	gmp_randinit_default(grt); //设置随机数生成算法为默认
	gmp_randseed_ui(grt, time(NULL)); //设置随机化种子为当前时间

	mpz_t key_p, key_q;
	mpz_init(key_p); //初始化q和p大素数
	mpz_init(key_q);

	mpz_urandomb(key_p, grt, sz);//随机生成一个0-2^1024的一个数
	mpz_urandomb(key_q, grt, sz);

	size_t sz = mpz_sizeinbase(key_p, 2);
	cout << "p的位数是" << sz << endl;

	sz = mpz_sizeinbase(key_q, 2);
	cout << "q的位数是" << sz << endl;

	//mpz_t key_pp, key_qq;
	mpz_nextprime(key_p, key_p);  //使用GMP自带的素数生成函数
	mpz_nextprime(key_q, key_q);
	gmp_printf("%Zd\n\n", key_p);
	gmp_printf("%Zd\n\n", key_q);

	//计算p*q的值，并存放在n中
	mpz_t key_n;
	mpz_init(key_n);

	mpz_mul(key_n, key_p, key_q); //计算p*q

	//计算(p-1)*(q-1)并将值放在key_f 中
	mpz_t key_f;
	mpz_init(key_f);
	mpz_sub_ui(key_p, key_p, 1);
	mpz_sub_ui(key_q, key_q, 1);
	mpz_mul(key_f, key_p, key_q);

	//找出符合要求的e，即公钥为(e,n) e通常取3,17和65537三个值，我们直接取e=65537
	mpz_t key_e;
	mpz_init_set_ui(key_e, 65537);

	//输出公钥（e,n)
	gmp_printf("%s (%Zd, %Zd)\n\n\n", "public key is:", key_n, key_e);

	//求e的逆元，即ed mod (f)=1;
	//用gmp自带的求数论逆元函数对其进行求解
	mpz_t key_d;
	mpz_init(key_d);
	mpz_invert(key_d, key_e, key_f);
	gmp_printf("%s (%Zd, %Zd)\n\n\n", "private key is:", key_d, key_e); //输出私钥


	//将明文m进行加密 C=m^e mod n
	mpz_t M, C;
	mpz_init(C);
	mpz_init_set_ui(M, 1234);
	//mpz_get_str(M, 16, M);
	mpz_powm(C, M, key_e, key_n);  //使用GMP中的模幂计算函数 C=M^e mod n;
	gmp_printf("%s %Zd\n\n", "the cipertxt is", C);

	mpz_t zz, t1, t2;
	mpz_init(zz);
	mpz_init(t1);
	mpz_init(t2);
	int sta = clock();
	double c1, c2, c3;
	for (int i = 0; i < 100000; i++)
	{
		//mpz_powm_ui(zz, C, 2, key_n);
		mpz_urandomb(t1, grt, sz << 1);
		mpz_urandomb(t2, grt, sz << 1);
		mpz_mul(zz,t1,t2);
		mpz_mod(zz, zz, key_n);
	}
	c1 = (double)(clock() - sta) / CLOCKS_PER_SEC;
	printf("%lf\n", c1);

	sta = clock();
	for (int i = 0; i < 100000; i++)
	{
		mpz_urandomb(zz, grt, sz << 1);
		mpz_jacobi(zz, key_n);
	}
	c2 = (double)(clock() - sta) / CLOCKS_PER_SEC;
	printf("%lf\n", c2);

	sta = clock();
	for (int i = 0; i < 100000; i++)
	{
		mpz_urandomb(zz, grt, sz << 1);
		mpz_invert(zz, zz, key_n);
	}
	c3 = (double)(clock() - sta) / CLOCKS_PER_SEC;
	printf("%lf\n", c3);
	cout << 2 * 256 * (c2 * 2 + c3)/100000 << endl;//Cocks
	cout << 2 * 256 * (c1 * 2 + c3)/100000 << endl;//proposal
	cout << 2 * 256 * (c1 * 4)/100000 << endl;//zhao et al
	//解密函数算法 M=C^d mod n
	mpz_t M2;
	//mpz_init_set_str(C2,M2,16);
	mpz_init(M2);
	mpz_powm(M2, C, key_d, key_n);   //使用GMP中的模幂计算函数 M=C^d mod n;
	gmp_printf("%s %Zd\n\n", "the M2 is", M2);

	mpz_clear(key_q);
	mpz_clear(key_p);
	mpz_clear(M);
	mpz_clear(C);
	mpz_clear(M2);
	mpz_clear(key_n);
	mpz_clear(key_f);
	mpz_clear(key_d);
	return 0;
}*/
/*char Nstr[1000] = "0xcb5645c59c402b0edcf96cbd6a7308b64aac2f37a3c6f96be7c421c4b7f0a4adbdecd88cbea1128352fb21baae583fe4ceb3fc93c4905803ad3e9214ada050d5c0ff785a13a5c9157c3154ad8d7015a2d239fe13ef836d3279c5cd5dc96013ac40f372a9c9226d2f5fe73f312c56e11d9cdfbf9fb0db627ac1a752f5f0bd2b29";
char d0str[1000] = "0x4f77b72b04e6fb2d02e5a43edef4784a2e22df0d42bfc7c9093a58ec35eb21a11962103be960b0088d0cc2e0dfb473bc2ba0a22cea1c73997442c8fab5e4bad22cd131055b0382eb9264ad40ec8257abaff11b33b173ffd0168039bf40dc203eb325d884d2845fd2b5a37f41a0f64183db0c256c244500000000000000000000";
const int t = 80;
mpz_t one;
int main()
{
	freopen("r.txt", "w", stdout);
	mpz_init(one);
	mpz_init_set_ui(one, 1);
	//mpz_t power;
	//mpz_set(power, one);
	//mpz_init(power);
	//for (int i = 0; i < t; i++)
	//{
	//	mpz_mul_ui(power, power, 2);//2^t
	//}
	//生成大素数算法
	gmp_randstate_t grt;
	gmp_randinit_default(grt); //设置随机数生成算法为默认
	gmp_randseed_ui(grt, time(NULL)); //设置随机化种子为当前时间

	mpz_t key_n, key_p, key_q, n4;//4*n
	mpz_init(key_n);
	mpz_init(key_p);
	mpz_init(key_q);
	mpz_init(n4);

	mpz_init_set_str(key_n, Nstr, 0);
	mpz_mul_ui(n4, key_n, 4);
	mpz_t key_e;
	mpz_init_set_ui(key_e, 65537);

	mpz_t d0;
	mpz_init(d0);
	mpz_init_set_str(d0, d0str, 0);

	//输出公钥（e,n)
	gmp_printf("%s (%Zd, %Zd)\n\n\n", "public key is:", key_n, key_e);

	mpz_t k1, k2;
	mpz_init(k1);
	mpz_init(k2);
	mpz_init_set_ui(k1, 25612);
	mpz_init_set_ui(k2, 25613);

	mpz_t k, d, d1, d2, tmp, tmpd, r, phi, tt, ttsqrt, tttmp, mul;
	mpz_init(k);
	mpz_init(d);
	mpz_init(tmpd);
	mpz_init(d1);
	mpz_init(d2);
	mpz_init(tmp);
	mpz_init(r);
	mpz_init(phi);
	mpz_init(tt);
	mpz_init(ttsqrt);
	mpz_init(tttmp);
	mpz_init(mul);
	mpz_set(k, k1);
	for (; mpz_cmp(k, k2) <= 0; mpz_add(k, k, one))
	{
		mpz_invert(d1, key_e, k);

		mpz_sub(tmp, d0, d1);
		mpz_fdiv_q(tmp, tmp, k);//tmp=(d0-d1)/k

		long long range = (((long long)1 << t) / mpz_get_ui(k)) + 1;
		cout << range << endl;
		for (long long j = 0; j < range; j++)
		{
			if (j % 10000000 == 0)
			{
				cout << j << endl;
			}
			mpz_add_ui(d2, tmp, j);
			mpz_mul(d2, d2, k);
			mpz_add(d, d2, d1);
			mpz_set(tmpd, d);
			mpz_mul(d, d, key_e);
			mpz_sub_ui(d, d, 1);//ed-1
			mpz_mod(r, d, k);
			if (mpz_cmp_ui(r, 0) == 0)
			{
				mpz_divexact(phi, d, k);//d是ed-1
				mpz_add_ui(tt, key_n, 1);
				mpz_sub(tt, tt, phi);
				mpz_set(tttmp, tt);//tttmp是p+q
				mpz_mul(tt, tttmp, tttmp);
				mpz_sub(tt, tt, n4);
				mpz_sqrt(ttsqrt, tt);
				mpz_add(key_q, tttmp, ttsqrt);//2q
				mpz_sub(key_p, tttmp, ttsqrt);//2p
				mpz_mul(mul, key_p, key_q);
				if (mpz_cmp(mul, n4) == 0)
				{
					mpz_div_ui(key_p, key_p, 2);
					mpz_div_ui(key_q, key_q, 2);
					gmp_printf("%Zd\n", key_p);
					gmp_printf("%Zd\n", key_q);
					gmp_printf("%Zd\n", tmpd);
					break;
				}
			}
		}
	}

	mpz_clear(key_q);
	mpz_clear(key_p);
	mpz_clear(key_e);
	mpz_clear(n4);
	mpz_clear(k);
	mpz_clear(key_n);
	mpz_clear(d);
	mpz_clear(d1);
	mpz_clear(d2);
	mpz_clear(tmp);
	mpz_clear(tmpd);
	mpz_clear(r);
	mpz_clear(phi);
	mpz_clear(tt);
	mpz_clear(ttsqrt);
	mpz_clear(tttmp);
	mpz_clear(mul);
	mpz_clear(k1);
	mpz_clear(k2);
	mpz_clear(d0);
	mpz_clear(one);
	return 0;
}
*/
char Nstr[1000] = "0xd46dd141810786e451320ca452b379024fd263501ae767760f3dcf34b79806b85e36b0fee538dac61a5872c37d051a8a026384d09f12b7e1adae7eb15c4d75878007ee0043c2186cf8999c59eb66f689f55baf190bd80e70bf47b553be76bd4efffc782a51b43314d54b83fc19461e1beb6021164f64723b505e5a619cb62335";
char estr[1000] = "0x92fbeeef2d40eb125234cfe4c063c4607f12aec7e3014b32fb4600e58c4eac1ec485192a1b03745632f2966311ad68bd1e49dd9d08b2bff67f58e214c8d7bae0142559994c24e347ff7555c86aa30ccd03cf794e6f00eead7f15e24f33da61fae11ec81e4e09bcc76c1a0ed5ca8c2f512856cdb42470beee7111a2410188697d";
mpz_t tp[1000];
const int t = 410;
mpz_t one;
int main()
{
	//freopen("r.txt", "w", stdout);
	mpz_init(one);
	mpz_init_set_ui(one, 1);
	mpz_t power;
	mpz_set(power, one);
	mpz_init(power);
	for (int i = 0; i < t; i++)
	{
		mpz_init(tp[i]);
		mpz_set(tp[i], power);
		mpz_mul_ui(power, power, 2);//2^t
	}
	//生成大素数算法
	gmp_randstate_t grt;
	gmp_randinit_default(grt); //设置随机数生成算法为默认
	gmp_randseed_ui(grt, time(NULL)); //设置随机化种子为当前时间

	mpz_t key_n, key_p, key_q, n4;//4*n
	mpz_init(key_n);
	mpz_init(key_p);
	mpz_init(key_q);
	mpz_init(n4);

	mpz_init_set_str(key_n, Nstr, 0);
	mpz_mul_ui(n4, key_n, 4);

	mpz_t key_e, d0;
	mpz_init_set_str(key_e, estr, 0);
	mpz_init(d0);
	//输出公钥（e,n)
	gmp_printf("%s (%Zd, %Zd)\n\n\n", "public key is:", key_n, key_e);

	mpz_t k1, k2;
	mpz_init(k1);
	mpz_init(k2);
	mpz_t k, d, d1, d2, tmp, tmpd, r, phi, tt, ttsqrt, tttmp, mul;
	mpz_init(k);
	mpz_init(d);
	mpz_init(tmpd);
	mpz_init(d1);
	mpz_init(d2);
	mpz_init(tmp);
	mpz_init(r);
	mpz_init(phi);
	mpz_init(tt);
	mpz_init(ttsqrt);
	mpz_init(tttmp);
	mpz_init(mul);
	mpz_set(k, k1);
	for (int i = 91; i < 400; i++)
	{
		for (int j = i + 1; j < 400; j++)
		{
			for (int k = j + 1; k < 400; k++)
			{
				for (int l = k + 1; l < 400; l++)
				{
					mpz_set(d0, tp[399]);
					mpz_add(d0, d0, tp[i - 1]);
					mpz_add(d0, d0, tp[j - 1]);
					mpz_add(d0, d0, tp[k - 1]);
					mpz_add(d0, d0, tp[l - 1]);

					mpz_set(tmpd, d);
					mpz_mul(d, d, key_e);
					mpz_sub_ui(d, d, 1);//ed-1
					mpz_mod(r, d, k);
					if (mpz_cmp_ui(r, 0) == 0)
					{
						mpz_divexact(phi, d, k);//d是ed-1
						mpz_add_ui(tt, key_n, 1);
						mpz_sub(tt, tt, phi);
						mpz_set(tttmp, tt);//tttmp是p+q
						mpz_mul(tt, tttmp, tttmp);
						mpz_sub(tt, tt, n4);
						mpz_sqrt(ttsqrt, tt);
						mpz_add(key_q, tttmp, ttsqrt);//2q
						mpz_sub(key_p, tttmp, ttsqrt);//2p
						mpz_mul(mul, key_p, key_q);
						if (mpz_cmp(mul, n4) == 0)
						{
							mpz_div_ui(key_p, key_p, 2);
							mpz_div_ui(key_q, key_q, 2);
							gmp_printf("%Zd\n", key_p);
							gmp_printf("%Zd\n", key_q);
							gmp_printf("%Zd\n", tmpd);
							break;
						}
					}
				}
			}
		}
	}
	for (; mpz_cmp(k, k2) <= 0; mpz_add(k, k, one))
	{
		mpz_invert(d1, key_e, k);

		mpz_sub(tmp, d0, d1);
		mpz_fdiv_q(tmp, tmp, k);//tmp=(d0-d1)/k

		long long range = (((long long)1 << t) / mpz_get_ui(k)) + 1;
		cout << range << endl;
		for (long long j = 0; j < range; j++)
		{
			if (j % 10000000 == 0)
			{
				cout << j << endl;
			}
			mpz_add_ui(d2, tmp, j);
			mpz_mul(d2, d2, k);
			mpz_add(d, d2, d1);
			mpz_set(tmpd, d);
			mpz_mul(d, d, key_e);
			mpz_sub_ui(d, d, 1);//ed-1
			mpz_mod(r, d, k);
			if (mpz_cmp_ui(r, 0) == 0)
			{
				mpz_divexact(phi, d, k);//d是ed-1
				mpz_add_ui(tt, key_n, 1);
				mpz_sub(tt, tt, phi);
				mpz_set(tttmp, tt);//tttmp是p+q
				mpz_mul(tt, tttmp, tttmp);
				mpz_sub(tt, tt, n4);
				mpz_sqrt(ttsqrt, tt);
				mpz_add(key_q, tttmp, ttsqrt);//2q
				mpz_sub(key_p, tttmp, ttsqrt);//2p
				mpz_mul(mul, key_p, key_q);
				if (mpz_cmp(mul, n4) == 0)
				{
					mpz_div_ui(key_p, key_p, 2);
					mpz_div_ui(key_q, key_q, 2);
					gmp_printf("%Zd\n", key_p);
					gmp_printf("%Zd\n", key_q);
					gmp_printf("%Zd\n", tmpd);
					break;
				}
			}
		}
	}

	mpz_clear(key_q);
	mpz_clear(key_p);
	mpz_clear(key_e);
	mpz_clear(n4);
	mpz_clear(k);
	mpz_clear(key_n);
	mpz_clear(d);
	mpz_clear(d1);
	mpz_clear(d2);
	mpz_clear(tmp);
	mpz_clear(tmpd);
	mpz_clear(r);
	mpz_clear(phi);
	mpz_clear(tt);
	mpz_clear(ttsqrt);
	mpz_clear(tttmp);
	mpz_clear(mul);
	mpz_clear(k1);
	mpz_clear(k2);
	mpz_clear(d0);
	mpz_clear(one);
	for (int i = 0; i < t; i++)
	{
		mpz_clear(tp[i]);
	}
	return 0;
}
// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
