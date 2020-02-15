// gnuplot_generate.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//这个函数的作用是生成全部每个画图文件的shell命令，到一个txt文件中然后再用gnuplot去调用

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cctype>
#pragma warning(disable : 4996)
char a[120][3] = { " ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
using namespace std;
const double pi = 3.14159;
const double e = 2.71828;
double gaosi(double x, double junzhi, double fangcha);
int main()
{
	
	int i = 0, j = 0, k = 0;
	//首先先针对得到的信息，产生数据文件，格式两列对应x和y，同时避免产生无效文件
	string data_path = "C:\\Users\\王志\\Desktop\\graph_series\\veigena_cell\\gaosi_data\\";
	string shell[20] = { "" };
	//开始自定义命令
	shell[0] = "set terminal postscript eps color solid  lw 2  \"Times_New_Roman\" 20";
	shell[1] = "set title ";//注意这里需要跟上title，可变部分
	shell[2] = "set xlabel \"Distances\"";
	shell[3] = "set ylabel \"Gauss sum\"";
	shell[4] = "set xrange [0.0:7.0]";
	shell[5] = "set mxtics 2";
	shell[6] = "set xtics 0, 0.25, 7";
	shell[7] = "set output "; //这里需要注意，加上输出的文件名
	shell[8] = "plot ";//后面加上数据文件名
	shell[9] = " using 1:2 with points linecolor -1 pointtype 7 pointsize 1,\"\" using 1:2 smooth csplines linecolor 7 linewidth 3";
	shell[10] = "set output ";

	shell[11] = "set term wxt";
	shell[12] = "set size ratio 0.8";
	shell[13] = "unset key";
	shell[14] = "set yrange [0.0:100]";

	FILE* in = fopen("get_peiduiname.txt", "r");
	if (in == NULL)
	{
		cout << "i can not find the file!" << endl;
		cin.get();
	}
	
	FILE* out = fopen("shell_all.gnuplot", "w");
	while (!feof(in))
	{
		char file_name[20] = {'\0'};
		fscanf(in, "%s\n", &file_name);
		int chang = strlen(file_name);
		for (int i = 0; i < 4; i++)
		{
			file_name[chang - i-1]='\0';

		}	
		//cout << "the file_name is:" << file_name << endl;
		//cin.get();
		fprintf(out, "%s\n", shell[0].c_str());
		fprintf(out, "%s\n", shell[12].c_str());
		fprintf(out, "%s\n", shell[13].c_str());
		//设置文件名		
		fprintf(out, "%s\n", (shell[1] + "\"" + file_name + "\"").c_str());
		//cout << "now the shell for 1 is:" << shell[1] + "\"" + file_name + "\"" << endl;
		//cin.get();
		//设置距离，刻度，坐标轴名字等等
		fprintf(out, "%s\n", shell[2].c_str());
		fprintf(out, "%s\n", shell[3].c_str());
		fprintf(out, "%s\n", shell[4].c_str());
		fprintf(out, "%s\n", shell[14].c_str());
		fprintf(out, "%s\n", shell[5].c_str());
		fprintf(out, "%s\n", shell[6].c_str());
		//设置输出文件名
		string out_name;
		out_name = shell[7] + "\"" + file_name + ".eps\"";
		fprintf(out, "%s\n", out_name.c_str());
		fprintf(out, "%s\n", (shell[8]+"\""+file_name+".txt\""+shell[9]).c_str());
		//cout << "now the shell for 2 is:" << shell[8] + "\"" + file_name + ".txt\"" + shell[9] << endl;
		//cin.get();
		fprintf(out, "%s\n", shell[10].c_str());
		fprintf(out, "%s\n", shell[11].c_str());
		fprintf(out, "\n");
	
	}
	fclose(out);
	fclose(in);


	cout << "all total work done!\n"; 
	cin.get();
	return 0;
}


double gaosi(double x, double junzhi, double fangcha)
{

	return (1 / (fangcha * pow(2 * pi, 0.5))) * pow(e, (-(x - junzhi) * (x - junzhi)) / (2 * pow(fangcha, 2)));
}