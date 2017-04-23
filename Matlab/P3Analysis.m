%{
The Pearson Type III Distribution Analysis Script
Matlab Script
Chinese/GB18030

MIT License

Copyright (c) 2017 llz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

disp('The Pearson Type III Distribution Analysis Script');
disp('Copyright (c) 2017 梁岚峥');
fprintf('\n');
clear;


%% 全局参数
%{
	Global_ImportData	导入经验数据
	Global_UseDefaultDataFile	使用默认数据文件
%}
Global_ImportData = true;	% <全局>导入经验数据
Global_UseDefaultDataFile = true;	% <全局>使用默认数据文件


%% 读取数据文件
%{
	DataFile	数据文件(标识号)
	DataFileName	数据文件名
	DataCount	数据组数
	DataRaw	原始数据
%}
if Global_ImportData
	DataFileName = 'data.txt';	% 数据文件名(默认值)
	if ~Global_UseDefaultDataFile
		disp('对话: 输入数据文件名:');
		DataFileName = input('$ ','s');	% 读取数据文件名(自定义)
	end
	DataFile = fopen(DataFileName,'r');	% 打开数据文件
	while DataFile < 0	% 打开数据文件错误处理
		fprintf('错误: 打开文件"%s"时发生了错误\n',DataFileName);
		disp('对话: 再次输入数据文件名:');
		DataFileName = input('$ ','s');
		DataFile = fopen(DataFileName,'r');
	end
	DataCount = fscanf(DataFile,'%d',[1]);	%读取数据组数
	DataRaw = fscanf(DataFile,'%f',[2,DataCount]);	%读取原始数据
	fclose(DataFile);	% 关闭数据文件
	fprintf('\n');
end


%% 绘制皮尔逊 III 型曲线
%{
	P3_Qbar	平均流量
	P3_Cv	变差系数
	P3_Cs	偏差系数
	P3_phip	离均系数
	P3_p	频率（x轴）
	P3_Qp	频率为p的洪峰流量（y轴）
	Delta_phip	计算误差时的离均系数
	Delta	误差
	Deltabar	平均误差
	DeltaCount	误差计数
	UserCommand	用户命令
%}
disp('对话: 输入平均流量 Qbar:');
UserCommand = input('$ ','s');	% 读取平均流量
P3_Qbar = str2double(UserCommand);
while isnan(P3_Qbar)	% 输入错误处理
	disp('错误: 输入应为实数');
	disp('对话: 再次输入平均流量 Qbar:');
	UserCommand = input('$ ','s');
	P3_Qbar = str2double(UserCommand);
end

disp('对话: 输入变差系数 Cv:');
UserCommand = input('$ ','s');	% 读取变差系数
P3_Cv = str2double(UserCommand);
while isnan(P3_Cv)	% 输入错误处理
	disp('错误: 输入应为实数');
	disp('对话: 再次输入变差系数 Cv:');
	UserCommand = input('$ ','s');
	P3_Cv = str2double(UserCommand);
end

disp('对话: 输入偏差系数 Cs:');
UserCommand = input('$ ','s');	% 读取偏差系数
P3_Cs = str2double(UserCommand);
while isnan(P3_Cs)	% 输入错误处理
	disp('错误: 输入应为实数');
	disp('对话: 再次输入偏差系数 Cs:');
	UserCommand = input('$ ','s');
	P3_Cs = str2double(UserCommand);
end

P3_p = linspace(0.001,0.999);	% 频率
if P3_Cs==0
	P3_phip = norminv(1-P3_p,0,1);
else
	P3_phip = 0.5*P3_Cs*gaminv(1-P3_p,4/P3_Cs^2,1)-2/P3_Cs;
end
P3_Qp = (1+P3_Cv*P3_phip)*P3_Qbar;
hold off;
plot(P3_p,P3_Qp);
hold on;
if Global_ImportData
	plot(DataRaw(2,:),DataRaw(1,:),'.');
	if P3_Cs==0
		Delta_phip = norminv(1-DataRaw(2,:),0,1);
	else
		Delta_phip = 0.5*P3_Cs*gaminv(1-DataRaw(2,:),4/P3_Cs^2,1)-2/P3_Cs;
	end
	Delta = (1+P3_Cv*Delta_phip)*P3_Qbar-DataRaw(1,:);
	DeltaCount = 0;
	Deltabar = 0;
	for i = 1:DataCount
		if ~(isnan(Delta(i))||isinf(Delta(i)))
			Deltabar = Deltabar + abs(Delta(i));
			DeltaCount = DeltaCount+1;
		end
	end
	Deltabar = Deltabar/DeltaCount;
	fprintf('结果: 绝对误差=%f\n',Deltabar);

	disp('对话: 根据经验点调整偏差系数？(yes/no)');
	UserCommand = input('$ ','s');	% 读取偏差系数
	while ~(strcmp(UserCommand,'yes')||strcmp(UserCommand,'no'))
		UserCommand = input('$ ','s');
	end

	if strcmp(UserCommand,'yes')
		MinDeltabar = inf;
		Bestm = 2;
		for m = 0:0.001:10
			P3_Cs = m*P3_Cv;
			if P3_Cs==0
				Delta_phip = norminv(1-DataRaw(2,:),0,1);
			else
				Delta_phip = 0.5*P3_Cs*gaminv(1-DataRaw(2,:),4/P3_Cs^2,1)-2/P3_Cs;
			end
			Delta = (1+P3_Cv*Delta_phip)*P3_Qbar-DataRaw(1,:);
			DeltaCount = 0;
			Deltabar = 0;
			for i = 1:DataCount
				if ~(isnan(Delta(i))||isinf(Delta(i)))
					Deltabar = Deltabar + abs(Delta(i));
					DeltaCount = DeltaCount+1;
				end
			end
			Deltabar = Deltabar/DeltaCount;
			if Deltabar<MinDeltabar
				MinDeltabar = Deltabar;
				Bestm = m;
			end
		end
		P3_Cs = Bestm*P3_Cv;
		if P3_Cs==0
			P3_phip = norminv(1-P3_p,0,1);
		else
			P3_phip = 0.5*P3_Cs*gaminv(1-P3_p,4/P3_Cs^2,1)-2/P3_Cs;
		end
		P3_Qp = (1+P3_Cv*P3_phip)*P3_Qbar;
		hold off;
		plot(DataRaw(2,:),DataRaw(1,:),'.');
		hold on;
		plot(P3_p,P3_Qp);
		fprintf('结果: 最优 m=%f Cs=%f 绝对误差=%f\n',Bestm,P3_Cs,MinDeltabar);
	end
end
