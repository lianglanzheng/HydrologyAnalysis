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
disp('Copyright (c) 2017 ����');
fprintf('\n');
clear;


%% ȫ�ֲ���
%{
	Global_ImportData	���뾭������
	Global_UseDefaultDataFile	ʹ��Ĭ�������ļ�
%}
Global_ImportData = true;	% <ȫ��>���뾭������
Global_UseDefaultDataFile = true;	% <ȫ��>ʹ��Ĭ�������ļ�


%% ��ȡ�����ļ�
%{
	DataFile	�����ļ�(��ʶ��)
	DataFileName	�����ļ���
	DataCount	��������
	DataRaw	ԭʼ����
%}
if Global_ImportData
	DataFileName = 'data.txt';	% �����ļ���(Ĭ��ֵ)
	if ~Global_UseDefaultDataFile
		disp('�Ի�: ���������ļ���:');
		DataFileName = input('$ ','s');	% ��ȡ�����ļ���(�Զ���)
	end
	DataFile = fopen(DataFileName,'r');	% �������ļ�
	while DataFile < 0	% �������ļ�������
		fprintf('����: ���ļ�"%s"ʱ�����˴���\n',DataFileName);
		disp('�Ի�: �ٴ����������ļ���:');
		DataFileName = input('$ ','s');
		DataFile = fopen(DataFileName,'r');
	end
	DataCount = fscanf(DataFile,'%d',[1]);	%��ȡ��������
	DataRaw = fscanf(DataFile,'%f',[2,DataCount]);	%��ȡԭʼ����
	fclose(DataFile);	% �ر������ļ�
	fprintf('\n');
end


%% ����Ƥ��ѷ III ������
%{
	P3_Qbar	ƽ������
	P3_Cv	���ϵ��
	P3_Cs	ƫ��ϵ��
	P3_phip	���ϵ��
	P3_p	Ƶ�ʣ�x�ᣩ
	P3_Qp	Ƶ��Ϊp�ĺ��������y�ᣩ
	Delta_phip	�������ʱ�����ϵ��
	Delta	���
	Deltabar	ƽ�����
	DeltaCount	������
	UserCommand	�û�����
%}
disp('�Ի�: ����ƽ������ Qbar:');
UserCommand = input('$ ','s');	% ��ȡƽ������
P3_Qbar = str2double(UserCommand);
while isnan(P3_Qbar)	% ���������
	disp('����: ����ӦΪʵ��');
	disp('�Ի�: �ٴ�����ƽ������ Qbar:');
	UserCommand = input('$ ','s');
	P3_Qbar = str2double(UserCommand);
end

disp('�Ի�: ������ϵ�� Cv:');
UserCommand = input('$ ','s');	% ��ȡ���ϵ��
P3_Cv = str2double(UserCommand);
while isnan(P3_Cv)	% ���������
	disp('����: ����ӦΪʵ��');
	disp('�Ի�: �ٴ�������ϵ�� Cv:');
	UserCommand = input('$ ','s');
	P3_Cv = str2double(UserCommand);
end

disp('�Ի�: ����ƫ��ϵ�� Cs:');
UserCommand = input('$ ','s');	% ��ȡƫ��ϵ��
P3_Cs = str2double(UserCommand);
while isnan(P3_Cs)	% ���������
	disp('����: ����ӦΪʵ��');
	disp('�Ի�: �ٴ�����ƫ��ϵ�� Cs:');
	UserCommand = input('$ ','s');
	P3_Cs = str2double(UserCommand);
end

P3_p = linspace(0.001,0.999);	% Ƶ��
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
	fprintf('���: �������=%f\n',Deltabar);

	disp('�Ի�: ���ݾ�������ƫ��ϵ����(yes/no)');
	UserCommand = input('$ ','s');	% ��ȡƫ��ϵ��
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
		fprintf('���: ���� m=%f Cs=%f �������=%f\n',Bestm,P3_Cs,MinDeltabar);
	end
end
