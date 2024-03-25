%Graph for Stray Field Comparison - Line Structures.

%Horizontal Lines

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 3);
% Specify sheet and range
opts.Sheet = "Sheet3";
opts.DataRange = "D26:F33";
% Specify column names and types
opts.VariableNames = ["LineThicknessum", "StrayFieldHAmleft", "UncertaintydeltaHAm"];
opts.VariableTypes = ["double", "double", "double"];
% Import the data
DielectricConstantVariablesS2 = readtable("C:\Users\prmee\OneDrive\Documents\UL 4th YR\Mar\06032024\Dielectric Constant Variables.xlsx", opts, "UseExcel", false);
% Clear temporary variables
clear opts



%Vertical Lines

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 3);
% Specify sheet and range
opts.Sheet = "Sheet3";
opts.DataRange = "D40:F45";
% Specify column names and types
opts.VariableNames = ["LineThicknessum", "StrayFieldHAmleft", "UncertaintydeltaHAm"];
opts.VariableTypes = ["double", "double", "double"];
% Import the data
DielectricConstantVariablesS1 = readtable("C:\Users\prmee\OneDrive\Documents\UL 4th YR\Mar\06032024\Dielectric Constant Variables.xlsx", opts, "UseExcel", false);
% Clear temporary variables
clear opts



plot(DielectricConstantVariablesS2,"LineThicknessum","StrayFieldHAmleft",LineStyle="--", Marker="X");
xlabel ('Line Thickness (\mum)','FontSize',10);
ylabel ('Magnitude of Stray Field, H (A/m)','FontSize',10)
title ('Magnitude of Stray Field with Changing Line Thickness','FontSize',10)
hold on
plot (DielectricConstantVariablesS1,"LineThicknessum","StrayFieldHAmleft",LineStyle="--", Marker="x")
hold off
legend('Horizontal Line Series','Vertical Line Series')



hxt = DielectricConstantVariablesS2(1:8,"LineThicknessum");
hx = table2array(hxt);
hyt = DielectricConstantVariablesS2(1:8,"StrayFieldHAmleft");
hy = table2array(hyt);
herrt = DielectricConstantVariablesS2(1:8,"UncertaintydeltaHAm");
herr = table2array(herrt);

vxt = DielectricConstantVariablesS1(1:6,"LineThicknessum");
vx = table2array(vxt);
vyt = DielectricConstantVariablesS1(1:6,"StrayFieldHAmleft");
vy = table2array(vyt);
verrt = DielectricConstantVariablesS1(1:6,"UncertaintydeltaHAm");
verr = table2array(verrt);



figure;
errorbar(hx,hy,herr,"vertical",Marker="x",LineStyle="--");
xlabel ('Line Thickness (\mum)','FontSize',10);
ylabel ('Magnitude of Stray Field, H (A/m)','FontSize',10)
hold on
e=errorbar(vx,vy,verr,"vertical",Marker="x",LineStyle="--");title ('Magnitude of Stray Field with Changing Line Thickness','FontSize',10)
xlabel ('Line Thickness (\mum)','FontSize',10);
ylabel ('Magnitude of Stray Field, H (A/m)','FontSize',10)
e.Color ='red';
hold off
legend('Horizontal Line Series','Vertical Line Series')