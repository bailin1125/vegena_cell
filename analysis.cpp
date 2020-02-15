//  vaerde_cell.cpp
//
//	input: the atom.config file
//	output: for specific atom in input file,you will get  the distance statistics files,by the vaerde cell method.

//
//  Created by 王志 on 2019/5/12.
//  Copyright  2019 王志. All rights reserved.
//
#include<iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <fstream>
#pragma warning(disable : 4996)
char style[20] = "style.ini";
using namespace std;
const double MAX = 7.0; //限制距离最大值，超过这个值无论是获取还是分析都会跳过						
void sort_plusxuhao(int *distance, int* xuhao, int num);
const double yuzhi = 0.2;
const double pi = 3.14159;
const double e = 2.71828;
const int yugu = 200;
double gaosi(double x, double junzhi, double fangcha);

char aa[120][3] = { " ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
double rule[114][114] = { 0.0 };
double dist[120][120] = { 0.0 };
int meatal_xuhao[89] = { 3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111 ,112 };

void get_style();
void amend_bs_rengong(void);
struct check
{
	int jishu;
	char elea[3];
	char eleb[3] ;
	double unsual;
	string file_name[yugu] ;
	double dis[yugu] ;
	bool flag ;//是否溢出的标志
};
//check结构用来批量储存异常结构
void fill_struct(struct check checka[]);
void generate_renwei_out(struct check checka[]);

int main(int argc, char*argv[])
{
	int i = 0, j = 0, k = 0;	
	int unusual_num = 0;
	get_style();
	
	//首先建立符合的三维数组来存储这些距离
	int ***save = new int **[114];
	for (i = 0; i < 114; i++)
	{
		save[i] = new int *[114];
		for (j = 0; j < 114; j++)
		{
			save[i][j] = new int[700];
		}
	}
	for (i = 0; i < 114; i++)
	{
		for (j = 0; j < 114; j++)
		{
			for (k = 0; k < 700; k++)
			{
				save[i][j][k] = 0;
			}
		}
	}
	char elementa[3];
	int a = 0;//序号a
	char elementb[3];
	int b = 0;//序号b
	double temp = 0.0;
	int location = 0;
	/*struct check checka[131];
	for (i = 0; i < 131; i++)
	{
		checka[i].flag = 0;
		checka[i].jishu = 0;
		checka[i].unsual = 0;
		for (j = 0; j < yugu; j++)
		{
			checka[i].dis[j] = 0;
			checka[i].file_name[j] = '\0';
		}

	}
	fill_struct(checka);*/

	//首先说明几个路径
	string path= "C:\\Users\\王志\\Desktop\\graph_series\\analysis_get_filename\\";//保存文件名的路径
	string path_aft = "_file_name.txt";
	string temp_path;
	string ori_file_path = "D:\\我的研究生\\学术\\相关项目\\weigena_cell\\vaerde_ori\\total\\";
	//string ori_file_path = "E:\\我的研究生\\学术\\相关项目\\vaerde_cell_ori\\total\\";
	string data_path = "C:\\Users\\王志\\Desktop\\graph_series\\veigena_cell\\data\\排除有机和合金\\";
	int jishu = 0;


	////不从原始数据获取再进行累加判断，直接读取目前的文件信息
	for (i = 1; i < 112; i++)
	{
		string xuhao;
		xuhao = aa[i];
		FILE* get_data = fopen((data_path + xuhao + ".txt").c_str(), "r");
		char temp[3] = { '\0' };
		//cout << "the find file_name is:" << (data_path + xuhao + ".txt").c_str() << endl;
		//cin.get();
		if (get_data == NULL)
		{
			cout << "i can not find the data file!" << endl;
			cin.get();
		}
		for (int j = i; j < 112; j++)
		{
			fscanf(get_data, "%s\t", &temp);
			//cout << "the temp is:" << temp<< endl;
			for (int k = 0; k < 700; k++)
			{
				fscanf(get_data, "%d\t", &save[i][j][k]);
			}
		}
		fclose(get_data);

	}
	//
	FILE* unsual = fopen("unsual_save", "wt");
	const char una[3] = { 'S' };
	const char unb[3] = { 'Cd' };
	const double unrule = 2.5;

	//int org_flag = 0;//对每个文件标记，是不是有机结构信息，是的话为1，并该文件直接跳过不再读取
	//int hejin_flag = 1;//每个文件标记，如果是合金的话就直接跳过
	//for (i = 1; i < 21; i++)
	//{
	//	cout << "now the :" << i << "file name save situation" << endl;
	//	temp_path.clear();
	//	char xuhao[20];
	//	sprintf(xuhao,"%d", i);
	//	temp_path = path + xuhao+path_aft;
	//	cout << "now the path is:" << temp_path.c_str() << endl;
	//	FILE* in = fopen(temp_path.c_str(), "r");		
	//	if (in == NULL)
	//	{
	//		cout << "i can not find the temp_path respectively,plesse" <<
	//			"consider " << endl;
	//		cout << "now the path is:" << temp_path.c_str() << endl;
	//		cin.get();
	//	}

	//	jishu = 0;
	//	while (!feof(in))//读取保存文件名字的文件
	//	{
	//		jishu++;
	//		if (jishu % 500 == 0)
	//			cout << "now the :" << jishu << "file\r" ;
	//		string file_name;			
	//		fscanf(in, "%s\n", file_name.c_str());
	//		char xuhao_a[20];
	//		string real_path = ori_file_path + xuhao + "\\"+file_name.c_str();
	//		FILE* get = fopen(real_path.c_str(), "r");
	//		//cout << "now the path is:" << real_path.c_str() << endl;
	//		if (get == NULL)
	//		{
	//			cout << "i can not find the file_name respectively,plesse" <<
	//				"consider " << endl;
	//			cout << "now the path is:" << real_path.c_str() << endl;
	//			cin.get();
	//		}

	//		int mm = 0;
	//		while (!feof(get))//在这个循环里确定是不是有机结构
	//		{
	//			char elementa[3];
	//			char elementb[3];
	//			fscanf(get, "%s\t", &elementa);
	//			//
	//			/*cout << file_name.c_str()<< endl;
	//			if (strcmp(file_name.c_str(),"111135")==0)
	//			{
	//				printf("%s\n", elementa);
	//			}*/
	//			fscanf(get, "%s\t", &elementb);
	//			//printf("%s\n", elementb);
	//			fscanf(get, "%lf\n", &temp);
	//			if ((strcmp(elementa, aa[1]) == 0 && strcmp(elementb, aa[6]) == 0) || (strcmp(elementa, aa[6]) == 0 && strcmp(elementb, aa[1]) == 0))
	//			{
	//				//cout << "got the organic structure!:" << file_name.c_str() << endl;
	//				//cin.get();
	//				org_flag = 1;
	//				break;
	//			}  
	//			for (mm = 0; mm < 89; mm++)
	//			{
	//				if (strcmp(elementa, aa[meatal_xuhao[mm]]) == 0)
	//					break;
	//			}
	//			if (mm == 89)
	//			{
	//				hejin_flag = 0;
	//				break;
	//			}

	//		}			
	//		fseek(get, 0, SEEK_SET);			
	//		/*if (hejin_flag == 1)
	//		{
	//			cout << "got meatal"<<file_name.c_str() << endl;
	//			cin.get();
	//		}*/
	//		if (org_flag == 1  || hejin_flag==1)
	//		{
	//			org_flag = 0;
	//			hejin_flag = 1;
	//			fclose(get);
	//			continue;
	//		}
	//         //在这里我们插入一个功能，要求记录上一些unsual的结构信息，存储下结构名字
	//		while (!feof(get))
	//		{
	//			int i = 0;				
	//			fscanf(get, "%s\t", &elementa);
	//			//printf("%s\n", elementa);
	//			fscanf(get, "%s\t", &elementb);
	//			//printf("%s\n", elementb);
	//			fscanf(get, "%lf\n", &temp);

	//			//这里先写好对于unsual的储存条件,注意这里是筛选某种或者某几种条件
	//			//if (strcmp(elementa, una) == 0 && strcmp(elementb, unb) == 0 && temp > unrule)
	//			//{
	//			//	fprintf(unsual, "%s\n", file_name.c_str());
	//			//	cout << "the file name is :" << file_name.c_str() << endl;
	//			//	cin.get();
	//			//	unusual_num++;
	//			//	if (unusual_num == 50)
	//			//	{
	//			//		cout << "has got 50 unusual file name!" << endl;
	//			//		
	//			//		cin.get();
	//			//		return 0;
	//			//	}
	//			//	break;
	//			//	//cin.get();
	//			//}

	//			//这里也是筛选特殊结构，不过是筛选批量化的
	//			//for (int i = 0; i < 131; i++)
	//			//{
	//			//	if (checka[i].flag==0 && strcmp(elementa, checka[i].elea) == 0 && strcmp(elementb, checka[i].eleb) == 0 && temp > checka[i].unsual)
	//			//	{						
	//			//		if ((checka[i].jishu+1) >= yugu)
	//			//		{
	//			//			checka[i].flag = 1;
	//			//			break;
	//			//		}
	//			//		else if (checka[i].file_name[checka[i].jishu-1] == file_name  && checka[i].jishu!=0)
	//			//		{
	//			//			break;
	//			//		}
	//			//		else
	//			//		{	
	//			//			checka[i].jishu++;
	//			//			checka[i].file_name[checka[i].jishu-1] = file_name.c_str();
	//			//			//cout << "the file name is:" << file_name.c_str() << endl;
	//			//			checka[i].dis[checka[i].jishu-1] = temp;
	//			//			//cout << "for " << elementa << " and" << elementb << " it's " << temp << endl;							
	//			//			//cin.get();
	//			//		}			

	//			//	}
	//			//}

	//			if (temp >= MAX)
	//			{
	//				//cout << "the distance is larger than :" << MAX << endl;
	//				//cout << "the temp is:" <<temp<< endl;
	//				//cin.get();
	//				continue;
	//			}
	//			//printf("%lf\n", temp);
	//			while (strcmp(elementa, aa[i]) != 0 && (i < 114))
	//			{
	//				i++;
	//			}
	//			if (i == 114)
	//			{
	//				//cout << "beyong the range!please check!" << endl;
	//				//cout << "the element a is:" << elementa << endl;
	//				//cout << "the file is:" << file_name.c_str() << endl;
	//				//cin.get();
	//				continue;
	//			}
	//			a = i;
	//			i = 0;
	//			while (strcmp(elementb, aa[i]) != 0 && (i < 114))
	//			{
	//				i++;
	//			}
	//			if (i == 114)
	//			{
	//				//cout << "beyong the range!please check!" << endl;
	//				//cout << "the elementb is:" << elementb << endl;
	//				//cout << "the file is:" << file_name.c_str() << endl;
	//				//cin.get();
	//				continue;
	//			}
	//			b = i;
	//			location = temp / 0.01;
	//			if (location > 699)
	//				location = 699;
	//			if (a < b)
	//			{
	//				save[a][b][location]++;
	//			}
	//			else
	//			{
	//				save[b][a][location]++;
	//			}
	//		}
	//		fclose(get);
	//	}
	//	fclose(in);
	//}
	//fclose(unsual);

	//先生成特殊的储存信息

	//然后开始生成144个文件来保存了	
	
	string hou=".txt";	
	for (i = 1; i < 114; i++)
	{
		string save_xuhao;
		char xuhao_save[20];
		strcpy(xuhao_save, aa[i]);
		save_xuhao = xuhao_save+hou;
		cout << "the save_file name is:" << save_xuhao.c_str() << endl;

		if (aa[i][0] == ' '|| aa[i][0] == '\0')
			break;	
		
		FILE* out;
		out = fopen((data_path+save_xuhao).c_str(), "wt");		
		for (j = i; j < 114; j++)
		{
			if (aa[j][0] == ' ' || aa[j][0] == '\0')
			{
				break;
			}				
			fprintf(out, "%s\t", aa[j]);			
			for (k = 0; k < 700; k++)
			{
				fprintf(out, "%d\t", save[i][j][k]);
			}
			fprintf(out, "\n");
		}
		fclose(out);
	}



	//现在我想简单生成每种组合的频率最多的五个距离及平均值
	//string pinlv_save = "result_every_3distances.txt";	
	//FILE* out = fopen(pinlv_save.c_str(), "wt");
	//int exit_flag = 0;
	//for (i = 1; i < 114; i++)
	//{
	//	if (aa[i][0] == ' ' || aa[i][0] == '\0')
	//		break;
	//	for (j = i; j < 114; j++)
	//	{
	//		if (aa[j][0] == ' ' || aa[j][0] == '\0')
	//		{				
	//			break;
	//		}
	//		int *cixu = new int[700];
	//		for (int k = 0; k < 700; k++)
	//		{
	//			cixu[k] = k;
	//		}
	//		sort_plusxuhao(save[i][j],cixu,700);

	//		double distance[10] = {0};
	//		double qiuhe = 0;
	//		int num = 1;
	//		for (int m = 0; m < 10; m++)
	//		{
	//			distance[m] = double(cixu[m] / 100.0);				
	//			//cout << distance[m] << endl;
	//		}
	//		qiuhe += distance[0];//这里我们默认频率最高的是最接近正确的
	//		for (int m = 1; m < 10; m++)
	//		{
	//			if ((distance[m]<distance[0] + yuzhi) && (distance[m] > distance[0] - yuzhi))
	//			{
	//				qiuhe += distance[m];
	//				num++;
	//			}
	//		}

	//		fprintf(out, "%s", aa[i]);
	//		fprintf(out, "-");
	//		fprintf(out, "%s\n", aa[j]);

	//		fprintf(out, "%lf\n", qiuhe/num);

	//		for (int m = 0; m < 5; m++)
	//		{
	//			fprintf(out, "%lf\t", distance[m]);
	//		}
	//		fprintf(out, "\n");
	//		delete[]cixu;
	//	}	
	//	
	//}
	//fclose(out);


	//generate_renwei_out(checka);
	//现在改变生成模式，建立多个文件，然后每个文件就是两列而已
	//k = 0;
	//string se_file_path = "C:\\Users\\王志\\Desktop\\graph_series\\veigena_cell\\独立元素配对data\\排除有机\\";
	//for (i = 1; i < 112; i++)
	//{
	//	for (j = i; j < 112; j++)
	//	{
	//		string se_filename;
	//		string befor,after;
	//		befor = aa[i];
	//		after = aa[j];
	//		se_filename = se_file_path+befor + "-" + after+ ".txt"; 
	//		while (save[i][j][k] == 0)
	//		{
	//			k++;
	//		}
	//		if (k == 700)
	//		{
	//			k = 0;
	//			continue;
	//		}
	//		FILE *out = fopen(se_filename.c_str(), "wt");
	//		for (k = 0; k < 700; k++)
	//		{
	//			fprintf(out, "%lf\t", double(k / 100.0));
	//			fprintf(out, "%d\n", save[i][j][k]);
	//		}
	//		fclose(out);
	//		k = 0;
	//	}
	//}

	//同时生成另外的一种文件，即高斯展开过后的数值
	k = 0;
	string gaosi_path = "C:\\Users\\王志\\Desktop\\graph_series\\veigena_cell\\gaosi_data\\50p_noorg_metal\\";
	for (i = 1; i < 112; i++)
	{
		for (j = i; j < 112; j++)
		{
			string gaosi_filename;
			string befor,after;
			befor = aa[i];
			after = aa[j];
			gaosi_filename = gaosi_path + befor + "-" + after + ".txt";
			while (save[i][j][k] == 0)
			{
				k++;
			}
			if (k == 700)
			{
				k = 0;
				//cout << "should jump this" << endl;
				continue;				
			}

			const int point_num = 50;//转换的点数
			const int jianxi_num = 8;//每两个点插多少个点
			const int chazhi_num = (point_num - 1)*jianxi_num + point_num;
			double temp[point_num] = { 0.0 };//每组数据我们转化为35个点
			double junzhi = 0.0;
			for (int m = 0; m < point_num; m++)
			{
				for (int n = 0; n < 700/point_num; n++)
				{
					junzhi += save[i][j][m * (700 / point_num) + n];
				}
				junzhi = junzhi / double(700 / point_num);
				temp[m] = junzhi;
				//在这里我们需要注意一个问题，为了结果好看我们将数据缩放到0-100之间
				
				/*for (int n = 0; n < 10; n++)
				{					
					temp[m] +=  10*gaosi(save[i][j][m * 10 + n], junzhi,1);
				}*/
			}
			double max = 0;
			double ratio = 1.0;
			for (int m = 0; m < point_num; m++)
			{
				if (max < temp[m])
					max = temp[m];
			}
			ratio = 100.0 / max;
			for (int m = 0; m < point_num; m++)
			{
				temp[m] = ratio * temp[m];
			}
			//这样得到了转化后的35个数值
			//然后进行插值
			double chazhi[chazhi_num] = { 0.0 };//插完后的点
			double cha[point_num - 1] = { 0.0 };//每个不同区间的增量
			for (int m = 0; m < (point_num-1); m++)
			{
				cha[m] = (temp[m + 1] - temp[m])/double(jianxi_num+1);
			}
			for (int m = 0; m < (point_num-1); m++)
			{
				for (int n = 0; n < (jianxi_num + 1); n++)
				{
					chazhi[m*(jianxi_num + 1) + n] = temp[m] + n * cha[m];
				}	

			}
			chazhi[chazhi_num-1] = temp[point_num - 1];

			//到这里获取了插值后的全部点，然后下面编写rule来获取截断点
			double up_num=0.0;
			double now_num = 0.0;
			const double yu = 0.038;
			const double max_ju = 3.6;//超过这个值认为不合理，需要削减
			const double temp_cer = 2.7;//如果按照rule没有定出距离值，则规定为这个值
			int rule_flag = 0;
			for (int m = 0; m < chazhi_num; m++)
			{
				if(chazhi[m]>15 && chazhi[m+1]<chazhi[m] && chazhi[m+5]<chazhi[m+4]&&chazhi[m+8]<chazhi[m+7])
				//找到了第一个下降峰
				{
					for (int n = m+1; n < chazhi_num; n++)
					{
						up_num = chazhi[n - 1];
						now_num = chazhi[n];
						if (now_num<75 &&now_num > (up_num-1.1) &&chazhi[n+4]>(chazhi[n+3]-0.8) && chazhi[n+6]>(chazhi[n+5]-0.8))
						{
							rule[i][j] = rule[j][i] = n*double(7.0/chazhi_num) + yu;
							while (rule[i][j] > max_ju)
							{
								rule[i][j] = rule[j][i] = rule[i][j] *0.9;
							}
							rule_flag = 1;
							break; 
						}
					}
				}
				if (rule_flag == 1)
					break;
			}
			if (rule[i][j]<1.0 || rule[i][j] == 0 || rule[i][j] == 0.0)
			{
				rule[i][j] = rule[j][i] = temp_cer;
			}
			

			FILE* out = fopen(gaosi_filename.c_str(), "w");
			if (out == NULL)
			{
				cout << "the file is null:" << gaosi_filename.c_str() << endl;
				cin.get();
			}
			for (int m = 0; m < chazhi_num; m++)
			{
				fprintf(out, "%lf\t", (7.0/chazhi_num)*m);
				fprintf(out, "%lf\n", chazhi[m]);
			}
			fclose(out);
			k = 0;
		}
	}




	//以上通过截断方式获得了机器规定的rule
	//下面需要基于人的规则进行更新或补充
	//1:对于H我们采用原始值
	for (int i = 1; i < 112; i++)
	{
		rule[1][i] = rule[i][1] = dist[1][i];
	}

	

	//3.对于缺省值我们规定为默认值1.6
	for (i = 0; i < 114; i++)
	{
		for (j = i; j < 114; j++)
		{
			if (rule[i][j] == 0.0)
			{
				rule[i][j] =rule[j][i]= 1.6;
			}
		}
	}
	//4.对于卤素Br35,I53进行检查，先用原始值去掉过大值
	const double max_ju = 3.5;//超过这个值认为不合理，需要削减
	for (i = 0; i < 114; i++)
	{
		if (rule[35][i] > max_ju)
		{
			rule[35][i] = rule[i][35] = rule[35][i] * 0.85;
		}
		if (rule[53][i] > max_ju)
		{
			rule[53][i] = rule[i][53] = rule[53][i] * 0.85;
		}
	}

	//然后最后检查一遍，看看还有没有超过限定值的距离
	//const double max_ju = 3.5;//超过这个值认为不合理，需要削减
	for (i = 0; i < 114; i++)
	{
		for (j = i; j < 114; j++)
		{
			while (rule[i][j] > max_ju)
			{
				rule[i][j] = rule[j][i] = rule[i][j] * 0.9;
			}
		}
	}


	//最后最后：对于人工筛选出来的，需要按照原始值或者使用新给定的值
	amend_bs_rengong();


	//cout <<"Ag I "<< rule[47][53] << endl;
	//cin.get();
	////开始生成判断的阈值文件
	ofstream rule_file;
	rule_file.open("rule_noorg_for_people_10-8.txt", ios::out);
	for (i = 1; i < 112; i++)
	{
		for (j = 1; j < 112; j++)
		{
			rule_file << aa[i] << '-' << aa[j] << ':';
			rule_file << rule[i][j] << '\t';
		}
		rule_file << endl;
	}
	rule_file.close();
	
	rule_file.open("rule_noorg_for_pc_10-8.txt", ios::out);
	for (i = 1; i < 112; i++)
	{
		for (j = 1; j < 112; j++)
		{
			rule_file << rule[i][j] << '\t';
		}
		rule_file << endl;
	}
	rule_file.close();

	//生成完了，开始删除空间
	for (i = 0; i < 114; i++)	
	{
		for (j = 0; j < 114; j++)
		{
			delete[] save[i][j];
		}		
	}
	delete[]save;
	cout << "all total work done!" << endl;
	cin.get();
	return 0;
   }


void sort_plusxuhao(int *distance, int* xuhao, int num)
{
	int i = 0, j = 0;
	double temp = 0.0;
	int temp_xuhao = 0;
	for (i = 0; i < num; i++)
	{
		for (j = i + 1; j < num; j++)
		{
			if (distance[i] < distance[j])
			{
				temp = distance[i];
				distance[i] = distance[j];
				distance[j] = temp;

				temp_xuhao = xuhao[i];
				xuhao[i] = xuhao[j];
				xuhao[j] = temp_xuhao;
			}
		}
	}
	return;
}
double gaosi(double x, double junzhi, double fangcha)
{

	return (1 / (fangcha * pow(2 * pi, 0.5))) * pow(e, (-(x - junzhi) * (x - junzhi)) / (2 * pow(fangcha, 2)));
}
void get_style()
{
	//char a[120][3] = { " ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	int i, j;
	char str[300];
	char qiana[3] = { '\0' };
	char houb[3] = { '\0' };
	char temp1[100];
	double rule = 0.0;
	int num1, num2;	
	int kk;
	FILE *in;
	in = fopen(style, "rt");
	if (in == NULL)
	{
		//cout<< "read the style.txt error!"<<endl;
		printf("error of reading style.txt!\n");
		return;
	}

	while (fgets(str, 300, in) != NULL)
	{
		if (strstr(str, "SBOND") != NULL)
		{
			break;
		}
	}
	for (i = 0; i < 986; i++)
	{
		fscanf(in, "%d", &kk);
		//printf("now is %d\n", kk);
		fscanf(in, "%s", qiana);
		fscanf(in, "%s", houb);
		fscanf(in, "%s", temp1);
		//fscanf(in, "%lf", rule);
		for (j = 0; j < 120; j++)
		{
			if (qiana[0] == aa[j][0] && qiana[1] == aa[j][1])
			{
				num1 = j;
				break;
			}
			//else
			//printf("can not find the qian\n");
		}
		for (j = 0; j < 120; j++)
		{
			if (houb[0] == aa[j][0] && houb[1] == aa[j][1])
			{
				num2 = j;
				break;
			}
		}
		fscanf(in, "%lf", &rule);
		//cout << "i got rule" << rule << endl;
		dist[num1][num2] = rule;
		dist[num2][num1] = dist[num1][num2];
		fgets(temp1, 100, in);

	}
	fclose(in);
	return;
}


void fill_struct(struct check checka[])
{
	int i = 0, j = 0;
	ifstream renwei;	
	renwei.open("rengong_unsual.txt", ios::in);
	if (!renwei.is_open())
	{
		cout << "i can not find the rengong_unsual file" << endl;
		cin.get();
	}

	for (i = 0; i < 131; i++)
		{
			 renwei>>checka[i].elea;
			 renwei >> checka[i].eleb;
			 renwei >> checka[i].unsual;
		}
	
	renwei.close();
	return;
}

void generate_renwei_out(struct check checka[])
{
	int i = 0, j = 0;
	ofstream fout;
	fout.open("check_filename.txt", ios::out);
	for (i = 0; i < 131; i++)
	{
		fout << checka[i].elea << '-' << checka[i].eleb << endl;
		for (j = 0; j < checka[i].jishu; j++)
		{
			fout << checka[i].file_name[j].c_str() << endl;
		}

	}
	
	fout.close();
	return;
}

void amend_bs_rengong(void)//通过此，把人找出的规则更新到最后的rule中
{
	int i = 0, j = 0;
	char elea[3], eleb[3];
	double temp = 0;

	ifstream fin1;//这个文件写明是需要用之前的数据
	fin1.open("choose_ori.txt", ios::in);
	if (!fin1.is_open())
	{
		cout << "i can not find the file for choose_ori!" << endl;
		cin.get();
	}
	
	while (fin1.good())
	{
		fin1 >> elea;
		fin1 >> eleb;
		for ( i = 0; i < 114; i++)
		{
			if (strcmp(elea, aa[i]) == 0)
			{
				break;
			}
		}
		if (i == 114)
		{
			break;
		}
		for (j = 0; j< 114; j++)
		{
			if (strcmp(eleb, aa[j]) == 0)
			{
				break;
			}
		}
		//cout << "elea is :" << elea;
		//cout << "the eleb is :" << eleb;
		//cin.get();
		rule[i][j] = rule[j][i] = dist[i][j];
	}
	fin1.close();

	ifstream fin2;//这个文件
	fin2.open("rengong_sure.txt", ios::in);
	if (!fin2.is_open())
	{
		cout << "i can not find the rengong_sure file!" << endl;
		cin.get();
	}
	while (fin2.good())
	{
		fin2 >> elea;
		fin2 >> eleb;
		fin2 >> temp;
		for (i = 0; i < 114; i++)
		{
			if (strcmp(elea, aa[i]) == 0)
			{
				break;
			}
		}
		if (i == 114)
			break;
		for (j = 0; j < 114; j++)
		{
			if (strcmp(eleb, aa[j]) == 0)
			{
				break;
			}
		}
		rule[i][j] = rule[j][i] = temp;
	}
	fin2.close();

	return;
}