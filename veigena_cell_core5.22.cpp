// vaerde_cell.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
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
#include <time.h>
#pragma warning(disable : 4996)
using namespace std;
//char plus_path[300];
//char wenjian[200];
//char outfile_path[300];
const bool generate_txt = false;
const double MAX = 7.0; //限制距离最大值，超过这个值无论是获取还是分析都会跳过						
const int cengshu = 3;							 //这个层数用来限制在相邻晶胞扩展时的层数
const int yanshen = cengshu * cengshu * cengshu; //说明扩了多少晶胞
char wenjian[200] = "C:\\Users\\DELL\\Desktop\\graph_series\\veigena_cell\\atom.config";
char outfile_path[300] = "C:\\Users\\DELL\\Desktop\\graph_series\\veigena_cell\\out1\\";
//char wenjian[200] = "/share/home/wangz/vaerde/atom.config";   
//char outfile_path[300] = "/share/home/wangz/vaerde/out1/";
char a[120][3] = { " ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
double dis(double *p1, double *p2)
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}
void sort(double *distance, int num);
class cell
{
public:
	cell(char *jiegou_name);
	int num;
	double **letice;
	double **p;
	double ***real_position;
	int *type;
};

cell::cell(char *name)
{
	int i, j, k;
	//cout << "expand the :" << cengshu << ":layer" << endl;
	char temp[300];
	double x_pian = 0.0;
	double y_pian = 0.0;
	double z_pian = 0.0;
	//strcpy(wenjian, "atom1.config");
	FILE *in;
	in = fopen(name, "rt");
	//system("pause");
	if (in == NULL)
	{
		printf("error of rading atom.config!\n");
		printf("the filename is :%s\n", wenjian);
		return;
	}
	fscanf(in, "%d", &num);

	type = (int *)malloc(num * sizeof(int));
	letice = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		letice[i] = (double *)malloc(3 * sizeof(double));
	}
	p = (double **)malloc(num * sizeof(double *));
	for (i = 0; i < num; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	real_position = (double ***)malloc(yanshen * sizeof(double **));
	for (i = 0; i < yanshen; i++)
	{
		real_position[i] = (double **)malloc(num * sizeof(double *));
		for (k = 0; k < num; k++)
			real_position[i][k] = (double *)malloc(3 * sizeof(double));
	}
	while (fgets(temp, 300, in) != NULL)
	{
		if (strstr(temp, "vector") != NULL)
			break;
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &letice[i][j]);
		}
	//x_pian = pow(pow(letice[0][0], 2) + pow(letice[0][1], 2) + pow(letice[0][2], 2), 0.5);
	//y_pian = pow(pow(letice[1][0], 2) + pow(letice[1][1], 2) + pow(letice[1][2], 2), 0.5);
	//z_pian = pow(pow(letice[2][0], 2) + pow(letice[2][1], 2) + pow(letice[2][2], 2), 0.5);
	fgets(temp, 300, in);
	fgets(temp, 300, in);
	for (i = 0; i < num; i++)
	{

		fscanf(in, "%d", &type[i]);
		fscanf(in, "%lf", &p[i][0]);
		fscanf(in, "%lf", &p[i][1]);
		fscanf(in, "%lf", &p[i][2]);
		fgets(temp, 300, in);
	}
	//int x_xishu = 0;
	//int y_xishu = 0;
	//int z_zishu = 0;

	for (i = 0; i < yanshen; i++)
	{

		for (j = 0; j < num; j++)
		{
			//x_xishu = i/3;
			//y_xishu = i / 3;
			//z_zishu = (i % 9) / 3;

			real_position[i][j][0] = letice[0][0] * p[j][0] + letice[1][0] * p[j][1] + letice[2][0] * p[j][2] + ((i % (cengshu * cengshu)) / cengshu - 1) * letice[0][0] + (i % cengshu - 1) * letice[1][0] - (i / (cengshu * cengshu) - 1) * letice[2][0];
			real_position[i][j][1] = letice[0][1] * p[j][0] + letice[1][1] * p[j][1] + letice[2][1] * p[j][2] + ((i % (cengshu * cengshu)) / cengshu - 1) * letice[0][1] + (i % cengshu - 1) * letice[1][1] - (i / (cengshu * cengshu) - 1) * letice[2][1];
			real_position[i][j][2] = letice[0][2] * p[j][0] + letice[1][2] * p[j][1] + letice[2][2] * p[j][2] + ((i % (cengshu * cengshu)) / cengshu - 1) * letice[0][2] + (i % cengshu - 1) * letice[1][2] - (i / (cengshu * cengshu) - 1) * letice[2][2];
		}
	}

	fclose(in);
}
class save
{
public:
	int point_num;
	int source_type; //记录从哪种元素找的最近邻
	//double **first_neighbour; //每个原子第一近邻的坐标
	int *type;		  //储存第一近邻原子的元素序号
	double *distance; //储存中心原子和第一近邻原子的距离
	//save(int num);
	save(void);
	~save(void);
};
save::save(void)
{
	int i, j, k;
	point_num = 0;
	source_type = 0;
	//num是第一近邻原子的个数
	//first_neighbour = new double *[num];
	/*for (i = 0; i < num; i++)
	{
		first_neighbour[i] = new double[3];
	}*/
	type = new int[30];
	distance = new double[30];
}//这个先默认存储了30个距离
save::~save(void)
{
	int i = 0, j = 0;
	delete[] distance;
	delete[] type;
	/*for (i = 0; i < num; i++)
	{
		delete[] first_neighbour[i];
	}
	delete[] first_neighbour;*/
	//cout << "the object is being destroied!" << endl;
}
int fullfill_face_function(double *face_function, double a[3], double b[3]);  //输入是两个点的坐标，输出结果是补充face_function
int san_face_topoint(double **face_function, double *point, double *in_face); //输入是三个平面的系数，以及第三个平面的一个点，输出是点的坐标或者判断平行之后false																	  //这个输入，第一个是中心原子的坐标，第二个是你计算出来的三面共点的坐标，第三个是便于找到的cell,以及找到的序号
int *judge_minum_distance(double *a, double *b, cell cell_a, int first, int two, int three, save *save_a, int *chose_point,int center);
void generate_outfile(char *name, save *save_a, cell cell_a);
void generate_txt_file(char* outfile_path); //用来生成7220个txt空文件
void fill_record_class(save save_a, cell cell_a, char *outfile_path);
void inv(int n, double **a, double **b);
double det(double **a, int n);
void generate_plus(save* save_a, cell cell_a, char* newfile_path);//根据相关结果生成增量文件，只说明增量信息
void sort_plusxuhao(double *distance, int* xuhao, int num);
int * judge_minum_distance_a(double *a, double *b, cell cell_a, int first, int two, int three, save *save_a, int *chose_point, int center); //������룬��һ��������ԭ�ӵ����꣬�ڶ������������������湲������꣬

int main(int argc, char *argv[])
{
	int i = 0, j = 0, k = 0, kk = 0, jj = 0, jjj = 0;
	int x1 = 0, x2 = 0;
	double temp = 0.0;
	double arrange = 0;
	double loop = 0;
	double save_infor = 0;
	double gene_file = 0;

	clock_t start;
	start = clock();
	clock_t finish;
	//从这里需要封装一层，变成可以批量处理的程序
	//首先我们先建立这么多txt
	//for (x1 = 0; x1 < argc; x1++)
	//{
	//	if (x1 == 1)
	//	{
	//		for (x2 = 0; x2 < strlen(argv[x1]); x2++)
	//		{
	//			wenjian[x2] = argv[x1][x2];
	//		}
	//		wenjian[strlen(argv[x1])] = '\0';
	//	}
	//	/*if (x1 == 2)
	//	{
	//		for (x2 = 0; x2 < strlen(argv[x1]); x2++)
	//		{
	//			outfile_path[x2] = argv[x1][x2];
	//		}
	//		outfile_path[strlen(argv[x1])] = '\0';
	//	}*/
	//	if (x1==2)
	//	{
	//		for (x2=0;x2<strlen(argv[i]);x2++)
	//		{
	//			plus_path[x2]=argv[x1][x2];
	//		}
	//		plus_path[strlen(argv[x1])]='\0';
	//	}


	//} 
	if (generate_txt == true)
	{
		//cout << "start to generate the out_file" << endl;
		generate_txt_file(outfile_path);
	}
	//cout << "generate 7220 txt file complete!" << endl;
	
	//char wenjian[200];
	//char outfile_path[300];
	
	//建立完成了这个用来储存的对象数组
	cell cell_a(wenjian);
	int num = cell_a.num * yanshen;
	double ** distance = new double*[num];
	for (j = 0; j < num; j++)
	{
		distance[j] = new double[num];
	}

	for (j = 0; j < num; j++)
	{
		for (k = j; k < num; k++)
		{
			distance[j][k] = dis(cell_a.real_position[j / cell_a.num][j%cell_a.num], cell_a.real_position[k / cell_a.num][k%cell_a.num]);
			distance[k][j] = distance[j][k];
		}
	}	
	save *save_a = new save[cell_a.num];
	//cout << "start the vaeede cell  processdure!" << endl;
	for (i = 0; i < cell_a.num; i++)
	{
		//cout << "start from the :" << i << "atom"
		//	 << "for total:" << num << endl;
		double **face_function = new double *[num];
		//建立num-1面的方程储存系数
		for (j = 0; j < num; j++)
		{
			face_function[j] = new double[4];
		}

		for (j = 0; j < num; j++) //任意其他原子与中心原子建立面方程
		{
			if (j != (((yanshen - 1) / 2)*cell_a.num+i))
			{
				if (fullfill_face_function(face_function[j], cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[j / cell_a.num][j % cell_a.num]) == 1)
				{
					//cout << "has written the " << j << "function" << endl;
				}
			}
			else
			{
				face_function[j][0] = face_function[j][1] = face_function[j][2] = face_function[j][3] = -10.0;
			}
		}
		//开始循环提取三个面的方程，来去输出得到的点的坐标（如果没有共点就不输出）
		//先建立一个存储三个面方程系数的东西，以及储存的标志0表示没有取出
		int *face_f_use = new int[num];
		for (j = 0; j < num; j++)
		{
			face_f_use[j] = 0;
		}

		for (x1 = 0; x1 < yanshen; x1++)
		{
			for (x2 = 0; x2 < i; x2++)
			{
				face_f_use[x1*cell_a.num + x2] = 2;
			}
		}
		//标记好对于之前找到的不再进行搜索
		double point[3];				//储存的共点坐标
		double in_face[3];				//储存三个面方程的d
		int chose_point[3];				//储存的是每个点找到的确定点的序号
		int *point_flag = new int[num]; //标记对中心原子来说，这个是不是最近邻原子,0表示不是
		for (j = 0; j < num; j++)
		{
			point_flag[j] = 0;
		}

		double **save_aa = new double*[3];
		for (j = 0; j < 3; j++)
		{
			save_aa[j] = new double[3];
		}
		double **save_bb = new double*[3];
		for (j = 0; j < 3; j++)
		{
			save_bb[j] = new double[3];
		}

		int ii = 0;
		int mm = 0;

		finish = clock();
		double arrange = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("the arrange memory has gone :%lf:seconds\n", arrange);		
		
		

		for (j = 0; j < num; j++) //每次提取的时候要把没用的那个避开
		{
			if (j != (((yanshen - 1) / 2)*cell_a.num + i) && face_f_use[j]!=2)
			{
				//printf("%d %d\n", i, j);
				face_f_use[j] = 1;
				temp = distance[(yanshen - 1) / 2 * cell_a.num + i][j];
				//temp = dis(cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[j / cell_a.num][j%cell_a.num]);
				if (temp > MAX)
				{
					face_f_use[j] = 0;
					continue;
				}
				for (k = j+1; k < num; k++)
				{
					if (k != (((yanshen - 1) / 2)*cell_a.num + i) && face_f_use[k] != 2)
					{
						face_f_use[k] = 1;
						temp = distance[(yanshen - 1) / 2 * cell_a.num + i][k];
						//temp = dis(cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[k / cell_a.num][k%cell_a.num]);
						if (temp > MAX)
						{
							face_f_use[k] = 0;
							continue;
						}
						for (kk = k+1; kk < num; kk++)
						{
							if (kk != (((yanshen - 1) / 2)*cell_a.num + i) && face_f_use[kk] != 2)
							{
								face_f_use[kk] = 1;								
								temp = distance[(yanshen - 1) / 2 * cell_a.num + i][kk];
								//temp = dis(cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[kk / cell_a.num][kk%cell_a.num]);
								if (temp > MAX)
								{
									face_f_use[kk] = 0;
									continue;
								}
								//遍历falg，如果提取出来了就放到我的save_face_function中
								for (jj = 0; jj < num; jj++)
								{
									if (face_f_use[jj] == 1)
									{
										//copy_double(save_face_funciton[jjj], face_function[jj]);
										for (x1 = 0; x1 < 3; x1++)
										{										
											//save_face_funciton[jjj][x1] = face_function[jj][x1];
											save_aa[jjj][x1] = face_function[jj][x1];											
										}
										in_face[jjj] = face_function[jj][3];
										jjj++;
										if (jjj == 3)
											break;
									}
								}
								/*if (j > 900)
								{
									printf("%d\t%d\t%d\t%d\n", i, j, k, kk);
									finish = clock();
									duration = (double)(finish - start) / CLOCKS_PER_SEC;
									printf("has gone:%lf,seconds\n", duration);
									cin.get();
									
								}	
								*/
								inv(3, save_aa, save_bb);
								point[0] = -save_bb[0][0] * in_face[0] - save_bb[0][1] * in_face[1] - save_bb[0][2] * in_face[2];
								point[1] = -save_bb[1][0] * in_face[0] - save_bb[1][1] * in_face[1] - save_bb[1][2] * in_face[2];
								point[2] = -save_bb[2][0] * in_face[0] - save_bb[2][1] * in_face[1] - save_bb[2][2] * in_face[2];
								if (point[0] == 0 && point[1] == 0 && point[2] == 0)
								{
									face_f_use[kk] = 0;
									jjj = 0;
									continue;
								}
								//printf("%lf %lf %lf\n", point[0], point[1], point[2]);
								/*
								if (san_face_topoint(save_face_funciton, in_face, point) == 1) ///这里需要注意没有解决in_face点的输入问题在，用的是第二个，输出是第三个
								{
									//cout << "i have succeeefuly save the point inforamtion" << endl;
									//cout << "random choose threee point's face_function" << j <<";"<< k <<";"<< kk << endl;
								}

								else //如果没有找到共点，那么就直接换一个面方程来试试
								{
									face_f_use[kk] = 0;
								jjj = 0;
									continue;
								}*/
								//printf("%lf %lf %lf\n", point[0], point[1], point[2]);
								
								int *pt = judge_minum_distance(cell_a.real_position[(yanshen - 1) / 2][i], point, cell_a, j, k, kk, save_a, chose_point,i);
								//int *pt = judge_minum_distance_a(cell_a.real_position[(yanshen - 1) / 2][i], point, cell_a, j, k, kk, save_a, chose_point, i); //������룬��һ��������ԭ�ӵ����꣬�ڶ������������������湲������꣬

								pt = chose_point;  
								////cout << *(pt) << endl;
								if (*(pt) != -100) //这个输入，第一个是中心原子的坐标，第二个是你计算出来的三面共点的坐标，第三个是便于找到的cell,输出是判断是否这个距离是最小的
								{
									//说明此时得到的点确实是合适的，此时应该将这个点对应的三个点储存到我们的save中
									// cout << "for " << i << "atom" << endl;
									point_flag[*(pt)] = point_flag[*(pt + 1)] = point_flag[*(pt + 2)] = 1;
									//cout << "three point are  valid:" << *(pt) << endl << *(pt + 1) << endl << *(pt + 2)   << endl;
									//cin.get();
								}
								face_f_use[kk] = 0;
								jjj = 0;
							}
						}
						face_f_use[k] = 0;
					}
				}
				face_f_use[j] = 0;
			}
		}
		//这个时候可以根据统计flag，来存储点的信息了
		finish = clock();
		loop = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("finish the loop has gone:%lf,seconds\n", loop);
		


		for (jj = 0; jj < num; jj++)
		{
			if (save_a[i].point_num == 30)
				break;
			if (point_flag[jj] == 1)
			{
				//temp = dis(cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[jj / cell_a.num][jj % cell_a.num]);
				temp = distance[(yanshen - 1) / 2 * cell_a.num + i][jj];
				if (temp < MAX)
					save_a[i].point_num++;				
			}
		}
		temp = 0.0;		
		save_a[i].source_type = cell_a.type[i];
		jjj = 0;
		for (jj = 0; jj < num; jj++)
		{
			if (point_flag[jj] == 1)		
			{				
				//temp = dis(cell_a.real_position[(yanshen - 1) / 2][i], cell_a.real_position[jj / cell_a.num][jj % cell_a.num]);
				temp = distance[(yanshen - 1) / 2 * cell_a.num + i][jj];
				if (temp > MAX)
					continue;				
				else
				{
					save_a[i].distance[jjj] = temp; 
					save_a[i].type[jjj] = cell_a.type[jj % cell_a.num];
					jjj++;
				}
			}
			if (jjj == 30)
				break;
		}

		delete[] point_flag;
		delete[] face_f_use;		
		for (j = 0; j < num; j++)
		{
			delete[] face_function[j];
		}
		delete[] face_function;		
		for (x1 = 0; x1 < 3; x1++)
		{
			delete[]save_aa[x1];
			delete[] save_bb[x1];
		}
		delete[]save_aa;
		delete[]save_bb;
		
		jjj = 0;

		finish = clock();
		save_infor = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("renew the loop has gone:%lf,seconds\n", save_infor);

	}	
	//从这里跳出来，说明完成了分析储存工作，开始要输出结果文件了 

	char plus_path[300] = "plus.txt";
	generate_outfile(plus_path, save_a, cell_a);
	//下面要根据每个获得的结果，进行结果文件的读写 	
	//for (i = 0; i < cell_a.num; i++)
	//{
		//fill_record_class(save_a[i], cell_a, outfile_path);
	//}		
	
	for (x1 = 0; x1 < num; x1++)
	{
		delete[] distance[x1];
	}
	delete[]distance;
	//system("pause");
	//delete []save_a;
	cout << "all total work done!" << endl;
	finish = clock();
	gene_file = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("save inforamtion and generate the file has gone:%lf,seconds\n", gene_file);	
	cin.get();
	return 0;
}

int * judge_minum_distance(double *a, double *b, cell cell_a, int first, int two, int three, save *save_a, int *chose_point,int center) //这个输入，第一个是中心原子的坐标，第二个是你计算出来的三面共点的坐标，
//第三个是便于找到的cell,输出是判断是否这个距离是最小的,(如果找到的距离最小值大于3个，视为有重复点的情况，则只返回和点的距离小的距离三个点)
//我们何不干脆
{
	//还需要加上一步，如果和这个距离相等的大于三个，说明有问题
	int i = 0, j = 0;
	int *pt = chose_point; //指针指向
	int num = cell_a.num * yanshen;
	double *distancea = new double[num];	
	int *xuhao = new int[num];
	int save_num = 0;

	for (i = 0; i < num; i++)
	{
		distancea[i] =  100;
		xuhao[i] = -1;
	}
	//初始化距离存储
	double temp = 0.0;
	double min = dis(a, b);
	//cout << "the distance is :" << min << endl;
	for (i = 0; i < num; i++)
	{
		temp = dis(b, cell_a.real_position[i / cell_a.num][i % cell_a.num]);			
		if (temp > MAX)
			continue;
		else if (temp < (min - 1e-4) && temp != 0.0)
		{			
			*(pt) = -100; //如果这个距离不是最小的，令第一个为-100，并返回
			//cout << chose_point[0] << endl;
			delete[]xuhao;
			delete[] distancea;			
			return pt;
		}
		else
		{
			distancea[j] = temp;
			xuhao[j] = i;
			j++;
		}	
		
	}

	num = j;
	//有j个距离了	
	sort_plusxuhao(distancea,xuhao, num);	
	j = 0;
	if (distancea[4] < min+1e-4  && distancea[4] > min - 1e-4)
	{
		int save_xuhao[9]; //先预存一个9个点的重复情况吧
		for (i = 0; i < 9; i++)
		{
			save_xuhao[i] = -1;
		}
		//cout << "has the same point situation!" << endl;
		//cin.get();
		//这里说明存在了重复点的情况 		
		for (i = 0; i < num; i++)
		{
			if (distancea[i] > min + 1e-3)
				break;
			else 
			{
				save_xuhao[j] = xuhao[i];
				j++;
				if (j == 9)
					break;
			}
		}		
		//获得了至多9个的相同距离的序号，下面判断找这些距离的最短三个,使用distance来储存距离并排序
		delete[] distancea;		
		delete[]xuhao;

		num = j;//记录了有几个共距离点
		double* dis_9 = new double[9];		
		for (i = 0; i < 9; i++)
		{
			 dis_9[i] = 100;
		}
		j = 0;
		for (i = 0; i < 9 ; i++)
		{
			if (save_xuhao[i] == -1)
				break;
			else
			{				
				dis_9[j]= dis(a, cell_a.real_position[save_xuhao[i] / cell_a.num][save_xuhao[i] % cell_a.num]);
				//dis_9[j] = distance[(yanshen - 1) / 2 * cell_a.num + center][save_xuhao[i]]; 
				j++;
			}
		}
		sort_plusxuhao(dis_9,save_xuhao,num);//排序之后我还是只要前三个
		*pt = save_xuhao[0];		
		*(pt + 1) = save_xuhao[1];
		*(pt + 2) = save_xuhao[2];
		/*int k = 0;
		for (i = 0; i < 3; i++)
		{
			for (k = 0; k < j; k++)
			{
				if (dis_9[i] == ori_9[k])
				{
					*(pt + i) = save_xuhao[k];
				}
			}
		}*/
		delete[] dis_9;		
		return pt;
	}
	else //这种情况是ok的
	{
		//这样的话我们应该得到了4个一样的最短距离
		//if ((distancea[0] < distancea[1]+1e-5) && (distancea[0] > distancea[1] - 1e-5)  &&(distancea[0] < distancea[2] + 1e-5) && (distancea[0] > distancea[2] - 1e-5))
		//{
		//	//cout << "yes i got it!" << endl;
		//	//cin.get(); 
		//}
		*(pt) = first;
		*(pt + 1) = two;
		*(pt + 2) = three;
		delete[] distancea;
		delete[]xuhao;
		return pt;
	}
}

int fullfill_face_function(double *face_function, double a[3], double b[3]) //输入是两个点的坐标，输出结果是补充face_function
{
	int i = 0, j = 0;
	int exit_flag = 0;
	double chuizhi[3];
	double in_face[3];
	for (i = 0; i < 3; i++)
	{
		chuizhi[i] = a[i] - b[i];
		in_face[i] = (a[i] + b[i]) / 2;
	}
	face_function[0] = chuizhi[0];
	face_function[1] = chuizhi[1];
	face_function[2] = chuizhi[2];
	face_function[3] = -(chuizhi[0] * in_face[0] + chuizhi[1] * in_face[1] + chuizhi[2] * in_face[2]);
	//cout << "the 4 paramters are :" << face_function[0]<<endl<<face_function[1] << endl <<face_function[2] << endl <<face_function[3] << endl;
	exit_flag = 1;
	return exit_flag;
}

int san_face_topoint(double **face_function, double *in_face, double *point) //输入是三个平面的系数，以及第三个平面的一个点，输出是点的坐标或者判断平行之后false
{
	//我怀疑没有传进来参数，看下吧
	int i, j;
	//for (i = 0; i < 3; i++)
	//{
	//	for (j = 0; j < 4; j++)
	//	{
	//		cout << "the paramter is :" << face_function[i][j] << endl;
	//	}
	//}
	int exit_flag = 0;
	//思路前两个面输出一个直线，再用直线和第三个面求交点
	double a = 0, b = 0, c = 0;
	double l = 0, m = 0, n = 0, o = 0, p = 0, q = 0, r = 0, s = 0, t = 0, c1 = 0, c2 = 0, x0 = 0, y0 = 0;
	l = face_function[0][0];
	m = face_function[0][1];
	n = face_function[0][2];
	o = face_function[1][0];
	p = face_function[1][1];
	q = face_function[1][2];
	r = face_function[2][0];
	s = face_function[2][1];
	t = face_function[2][2];
	c1 = face_function[0][3];
	c2 = face_function[1][3];
	double tolerance = 0.01;
	double Determinant = 0.0;
	Determinant = l * (p*t - q * s) - m * (o*t - q * r) + n * (o*s - p * r);
	//cout << "the hang lie shi velue is:" << Determinant << endl;

	if (Determinant<tolerance && Determinant>-tolerance)
	{
		//cout << "three face has parallel situation,so quit" << endl;
		return exit_flag;
	}
	else
	{
		//输出直线的方向向量
		a = m * q - n * p;
		b = n * o - l * q;
		c = l * p - m * o;
		//cout << "the line direction vector are:" << a << endl << b << endl << c << endl;
		//求出交线上面的一点，并且令z=0
		x0 = (p*c1 - m * c2) / (o*m - p * l);
		y0 = (l*((m*c2 - p * c1) / (o*m - p * l)) - c1) / m;
		//cout << "the x0 is:" << x0 << "and the y0 is:" << y0 << endl;
		//用行列式求解交点坐标，并且我们检查下这里面中的点取的对不对
		//cout << "the in face point is:" << in_face[0] << endl << in_face[1] << endl << in_face[2] << endl; 
		double z = 0.0;
		double fenmu = 0.0;
		fenmu = r * a + s * b + t * c;
		z = ((in_face[0] - x0) * r + (in_face[1] - y0) * s + (in_face[2]) * t) / fenmu;
		//cout << "check the fenmu is:" << r * a + s * b + t * c << endl;
		//cout << "the z is:" << z << endl;

		point[0] = x0 + a * z;
		point[1] = y0 + b * z;
		point[2] = c * z;
		if ((c<0.00001&&c>-0.00001) || (m<0.00001&&m>-0.00001)|| (fenmu<0.00001&&fenmu>-0.00001))//说明计算中出现了分母为0的情况
		{
			//if (x0 > 1e6||x0<-1e6)
				//cout << "too large" << endl;
			//point[0] = point[1] = point[2] = 100;
			return exit_flag;
		}
		exit_flag = 1;
		//cout << "the inface point is:" << point[0] << endl << point[1] << endl << point[2] << endl;
		//cin.get();
		return exit_flag;
	}
}

void generate_outfile(char *name, save *save_a, cell cell_a) //输入是你储存得到的信息(针对当前结构的输出），输出是在当前目录下的txt文件，记录结果用来检查
{
	int i = 0, j = 0;
	
	FILE *out = fopen(name, "wt");	
	for (i = 0; i < cell_a.num; i++)
	{
		
		for (j = 0; j < save_a[i].point_num; j++)
		{
			fprintf(out, "%s\t", a[save_a[i].source_type]);						
			fprintf(out, "%s\t", a[save_a[i].type[j]]);
			fprintf(out, "%lf\n", save_a[i].distance[j]);
		}
		//fprintf(out, "\n");
	}
	fclose(out);
	//cout << "has generated the file:" << name << endl;
	return;
}

void fill_record_class(save save_a, cell cell_a, char *outfile_path) //输入分别是save用来储存最近邻信息 ，结果是将数据填充到对应的txt中
{
	int i = 0, j = 0, k = 0;
	int save[700];
	char complete_path[300];	
	char file_name[50];
	char name_temp[20];
	int location; //根据距离确定他在哪个位置需要插入或者++
	int temp = 0;
	FILE *out;
	//这里需要注意，我们建立的txt名字需要对应好，小序号在前
	int before, after;	
	for (i = 0; i < save_a.point_num; i++)
	{
		strcpy(complete_path, outfile_path);
		if (save_a.source_type > save_a.type[i])
		{
			before = save_a.type[i];
			after = save_a.source_type;
		}
		else
		{
			before = save_a.source_type;
			after = save_a.type[i];
		}
		strcpy(file_name, a[before]);
		strcat(file_name, "-");
		strcpy(name_temp, a[after]);
		strcat(file_name, name_temp);
		strcat(complete_path, file_name);
		//cout << "i will write in the filename is :" << complete_path << endl;
		out = fopen(complete_path, "r");
		if (out == NULL)
		{
			//cout << "i can not find the file respectavely!" << endl;
			//cout << "the outfile_path is :" << complete_path<<endl;
			//cin.get();
		}
		if (save_a.distance[i] > MAX)
		{
			fclose(out);
			continue;
		}
		location = save_a.distance[i] / 0.01;
		for (j = 0; j < 700; j++)
		{
			fscanf(out, "%d\n", &save[j]);
		}
		save[location]++;
		fclose(out);
		out = fopen(complete_path, "wt");
		for (j = 0; j < 700; j++)
		{
			fprintf(out, "%d\n", save[j]);
		}
		fclose(out);
	}
	return;
}

void generate_txt_file(char* outfile_path)
{
	int i = 0, j = 0, k = 0,m=0;
	char path[300];
	strcpy(path, outfile_path);
	//cout << "the path is :" << path << endl;
	char file_name[300];
	char temp[20];
	FILE *out;
	for (i = 0; i < 119; i++)
	{
		for (j = i; j < 119; j++)
		{
			strcpy(file_name, a[i + 1]);
			strcpy(temp, a[j + 1]);
			strcat(file_name, "-");
			strcat(file_name, temp);
			strcat(path, file_name);
			//cout << "the path is :" << path << endl;
			out = fopen(path, "w");	
			//if (out == NULL)
			//{
			//	cout << "i can not find the path" << endl;
			//	cin.get();
			//}
			for (k = 0; k < 700; k++)
			{
				fprintf(out, "%d\n", m);
			}
			fclose(out);
			strcpy(path, outfile_path);
		}
	}
	//cout << "generate out_file finished" << endl;
	return;
}

void sort(double *distance, int num)
{
	int i = 0, j = 0;
	double temp = 0.0;
	for (i = 0; i < num; i++)
	{
		for (j = i + 1; j < num; j++)
		{
			if (distance[i] > distance[j])
			{
				temp = distance[i];
				distance[i] = distance[j];
				distance[j] = temp;
			}
		}
	}
	return;
}

void sort_plusxuhao(double *distance, int* xuhao,int num)
{
	int i = 0, j = 0;
	double temp = 0.0;
	int temp_xuhao = 0;
	for (i = 0; i < num; i++)
	{
		for (j = i + 1; j < num; j++)
		{
			if (distance[i] > distance[j])
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


double det(double **a, int n)
{
	int i, j, k;
	int ii, jj;
	int a_i, a_j;
	a_i = -1;
	a_j = -1;
	double m;
	m = 0;
	double **b;
	int f;
	b = (double **)malloc((n - 1) * sizeof(double*));
	for (i = 0; i < n - 1; i++)
	{
		b[i] = (double *)malloc((n - 1) * sizeof(double));
	}
	if (n == 1)
	{
		m = a[0][0];
	}
	else
	{
		j = 0;
		for (i = 0; i < n; i++)
		{
			f = 1;
			for (k = 0; k < i + j; k++)
			{
				f = f * (-1);
			}
			for (ii = 0; ii < n; ii++)
			{
				for (jj = 0; jj < n; jj++)
				{
					if (ii < i)
					{
						a_i = ii;
					}
					if (ii > i)
					{
						a_i = ii - 1;
					}
					if (jj < j)
					{
						a_j = jj;
					}
					if (jj > j)
					{
						a_j = jj - 1;
					}
					if (ii != i && jj != j)
					{
						b[a_i][a_j] = a[ii][jj];
					}
				}
			}
			m = m + a[i][j] * det(b, n - 1)*f;
		}
	}
	for (i = 0; i < n - 1; i++)
	{
		free(b[i]);
	}
	free(b);
	return m;
}
void inv(int n, double **a, double **b)
{
	int i, j, k;
	int ii, jj;
	int a_i, a_j;
	a_i = -1;
	a_j = -1;
	double det_a;
	double **c;
	int f;

	c = (double **)malloc((n - 1) * sizeof(double*));
	for (i = 0; i < n - 1; i++)
	{
		c[i] = (double *)malloc((n - 1) * sizeof(double));
	}

	det_a = det(a, n);
	if (fabs(det_a) < 0.000000001)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				b[i][j] = 0;
			}
		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				f = 1;
				for (k = 0; k < i + j; k++)
				{
					f = f * (-1);
				}
				for (ii = 0; ii < n; ii++)
				{
					for (jj = 0; jj < n; jj++)
					{
						if (ii < i)
						{
							a_i = ii;
						}
						if (ii > i)
						{
							a_i = ii - 1;
						}
						if (jj < j)
						{
							a_j = jj;
						}
						if (jj > j)
						{
							a_j = jj - 1;
						}
						if (ii != i && jj != j)
						{
							c[a_i][a_j] = a[ii][jj];
						}
					}
				}
				b[j][i] = det(c, n - 1) / det_a * f;
			}
		}
	}
	for (i = 0; i < n - 1; i++)
	{
		free(c[i]);
	}
	free(c);
	return;
}

void generate_plus(save* save_a, cell cell_a, char* newfile_path)
{
	int i = 0, j = 0;
	FILE* out;
	out = fopen(newfile_path, "wt");
	printf("the new_file path is:%s\n", newfile_path);
	for (i = 0; i < cell_a.num; i++)
	{
		printf("the point num is:%d\n", save_a[i].point_num);
		for (j = 0; j < save_a[i].point_num; j++)
		{
			//if (save_a[i].distance[j] < 1e-6)
				//continue;
			printf("the type are:%s\n", a[save_a[i].source_type]);
			fprintf(out, "%s-", a[save_a[i].source_type]);
			printf("ok?\n");

			printf("the type are:%d\n", save_a[i].type[j]);
			fprintf(out, "%s\t", a[save_a[i].type[j]]);
			printf("ok?\n");
			printf("the distance is:%lf\n", save_a[i].distance[j]);
			fprintf(out, "%lf\n", save_a[i].distance[j]);
			printf("ok?\n");
		}
	}
	fclose(out);
	return;
}
void copy_double(double *a, double *b) 
{int i = 0;
	for (i = 0; i < 4; i++)
	{
		a[i] = b[i];
	}
	return;
}




int * judge_minum_distance_a(double *a, double *b, cell cell_a, int first, int two, int three, save *save_a, int *chose_point, int center) //������룬��һ��������ԭ�ӵ����꣬�ڶ������������������湲������꣬
//�������Ǳ����ҵ���cell,������ж��Ƿ������������С��,(����ҵ��ľ�����Сֵ����3������Ϊ���ظ�����������ֻ���غ͵�ľ���С�ľ���������)
//���Ǻβ��ɴ�
{
	//����Ҫ����һ������������������ȵĴ���������˵��������
	int i = 0, j = 0;
	int *pt = chose_point; //ָ��ָ��
	int num = cell_a.num * yanshen;
	double *distance = new double[num];
	double *original = new double[num];
	int save_num = 0;

	double temp = 0.0;
	double min = dis(a, b);
	//cout << "the distance is :" << min << endl;
	for (i = 0; i < num; i++)
	{
		temp = dis(b, cell_a.real_position[i / cell_a.num][i % cell_a.num]);
		if (i == (((yanshen - 1) / 2)*cell_a.num + center))
		{
			temp = 100;
		}
		original[i] = temp;
		if (temp < (min - 0.001) && temp != 0.0)
		{
			//cout << "i find the liwai situation:" << temp << endl;
			*(pt) = -100; //���������벻����С�ģ����һ��Ϊ-100��������
			//cout << chose_point[0] << endl;
			delete[] distance;
			delete[] original;
			return pt;
		}
	}
	for (i = 0; i < num; i++)
	{
		distance[i] = original[i];
	}
	sort(distance, num);
	/*for (i = 0; i < num; i++)
	{
		cout << "distance after paixu is:" << endl;
		cout << distance[i] << endl;
	}
	cin.get();*/
	int save_xuhao[9]; //��Ԥ��һ��15������ظ������
	for (i = 0; i < 9; i++)
	{
		save_xuhao[i] = -1;
	}
	if (distance[3]< min+1e-4 && distance[3] > min - 1e-4)
	{
		//cout << "has the same point situation!" << endl;
		//cin.get();
		//����˵���������ظ������� 		
		for (i = 0; i < num; i++)
		{
			if (original[i] < (min + 0.0001) && original[i] > (min - 0.0001))
			{
				save_xuhao[j] = i;
				j++;
				if (j == 9)
					break;
			}
		}
		/*for (i = 0; i < num; i++)
		{
			cout << "distance after paixu is:" << endl;
			cout << distance[i] << endl;
		}
		cin.get();*/
		/*cout << "same point num is :" << j  << endl;
		for (i = 0; i < j; i++)
		{
			cout << distance[i] << endl;
		}
		cin.get();*/
		//���������15������ͬ�������ţ������ж�����Щ������������,ʹ��distance��������벢����
		j = 0;
		for (i = 0; i < 9; i++)
		{
			if (save_xuhao[i] == -1)
				break;
			else
			{
				original[j] = distance[j] = dis(a, cell_a.real_position[save_xuhao[i] / cell_a.num][save_xuhao[i] % cell_a.num]);
				j++;
			}
		}
		sort(distance, j - 1);//����֮���һ���ֻҪǰ����
		int k = 0;
		for (i = 0; i < 3; i++)
		{
			for (k = 0; k < j - 1; k++)
			{
				if (distance[i] == original[k])
				{
					*(pt + i) = save_xuhao[k];
				}
			}
		}
		delete[] distance;
		delete[] original;
		return pt;
	}
	else //���������ok��
	{
		//�����Ļ�����Ӧ�õõ���4��һ������̾���
		//if ((distance[0] < distance[1] + 1e-5) && (distance[0] > distance[1] - 1e-5) && (distance[0] < distance[2] + 1e-5) && (distance[0] > distance[2] - 1e-5))
		//{
		//	//cout << "yes i got it!" << endl;
		//	//cin.get(); 
		//}
		*(pt) = first;
		*(pt + 1) = two;
		*(pt + 2) = three;
		delete[] distance;
		delete[] original;
		return pt;
	}
}


	