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
    
    

