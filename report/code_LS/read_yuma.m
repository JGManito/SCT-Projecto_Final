function [dataset_out] = read_yuma(input_file)
%READ_YUMA Reads Yuma data and outputs a formated data set
%   This function receives as input a Yuma Almanac text file, obtained from
%   celestrak.com, and parses it.

%Open the file
fid=fopen(input_file);
%Check if file opened correctly
if fid == -1
    disp('Input file error');
    return;
end

line=fgetl(fid); %read the first line
i=0; %Initialize the line counter

while ischar(line)
    i=i+1;

    %Read the input file
    line=fgetl(fid);
    dataset_out(i).id=sscanf(line,"ID:                         %d\n");
    %disp(dataset_out(i).id);
    line=fgetl(fid);
    dataset_out(i).health=sscanf(line,"Health:                     %d\n");
    %disp(dataset_out(i).health);
    line=fgetl(fid);
    dataset_out(i).eccentricity=sscanf(line,"Eccentricity:               %f\n");
    %disp(dataset_out(i).eccentricity);
    line=fgetl(fid);
    dataset_out(i).time=sscanf(line,"Time of Applicability(s):  %f\n");
    %disp(dataset_out(i).time);
    line=fgetl(fid);
    dataset_out(i).inclination=sscanf(line,"Orbital Inclination(rad):   %f\n");
    %disp(dataset_out(i).inclination);
    line=fgetl(fid);
    dataset_out(i).omega_dot=sscanf(line,"Rate of Right Ascen(r/s):  %f\n");
    %disp(dataset_out(i).omega_dot);
    line=fgetl(fid);
    dataset_out(i).sqrt_A=sscanf(line,"SQRT(A)  (m 1/2):           %f\n");
    %disp(dataset_out(i).sqrt_A);
    line=fgetl(fid);
    dataset_out(i).omega_zero=sscanf(line,"Right Ascen at Week(rad):  %f\n");
    %disp(dataset_out(i).omega_zero);
    line=fgetl(fid);
    dataset_out(i).arg_perigee=sscanf(line,"Argument of Perigee(rad):   %f\n");
    %disp(dataset_out(i).arg_perigee);
    line=fgetl(fid);
    dataset_out(i).mean_anomaly=sscanf(line,"Mean Anom(rad):             %f\n");
    %disp(dataset_out(i).mean_anomaly);
    line=fgetl(fid);
    dataset_out(i).af0=sscanf(line,"Af0(s):                    %f\n");
    %disp(dataset_out(i).af0);
    line=fgetl(fid);
    dataset_out(i).af1=sscanf(line,"Af1(s/s):                  %f\n");
    %disp(dataset_out(i).af1);
    line=fgetl(fid);
    dataset_out(i).week=sscanf(line,"week:                        %d\n");
    %disp(dataset_out(i).week);
    
    %Skip the last line and the first line of the next block
    line=fgetl(fid);
    line=fgetl(fid);
    
end
    %---------DEBUG--------%
% for j=1:24
%     A = [272.847 268.126; 272.847 161.786; 272.847 11.676; 272.847 41.806; 332.847 80.956; 352.847 173.336; 332.847 309.976; 332.847 204.376; 32.847 111.876; 32.847 11.796; 32.847 339.666; 32.847 241.556; 92.847 135.226; 92.847 265.446; 92.847 35.156; 92.847 167.356; 152.847 197.046; 152.847 302.596; 152.847 333.686; 152.847 66.066; 212.847 238.886; 212.847 345.226; 212.847 105.206; 212.847 135.346];
%     dataset_out(j).omega_zero = A(j,1);
%     dataset_out(j).arg_perigee = A(j,2);
% end
% 
% dataset_out_ = dataset_out(:,1:24);
% dataset_out=[];
% dataset_out=dataset_out_;
    

