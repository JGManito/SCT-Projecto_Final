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
	
    line=fgetl(fid);
    dataset_out(i).health=sscanf(line,"Health:                     %d\n");
    
    line=fgetl(fid);
    dataset_out(i).eccentricity=sscanf(line,"Eccentricity:               %f\n");
    
    line=fgetl(fid);
    dataset_out(i).time=sscanf(line,"Time of Applicability(s):  %f\n");
    
    line=fgetl(fid);
    dataset_out(i).inclination=sscanf(line,"Orbital Inclination(rad):   %f\n");
    
    line=fgetl(fid);
    dataset_out(i).omega_dot=sscanf(line,"Rate of Right Ascen(r/s):  %f\n");
    
    line=fgetl(fid);
    dataset_out(i).sqrt_A=sscanf(line,"SQRT(A)  (m 1/2):           %f\n");
    
    line=fgetl(fid);
    dataset_out(i).omega_zero=sscanf(line,"Right Ascen at Week(rad):  %f\n");
    
    line=fgetl(fid);
    dataset_out(i).arg_perigee=sscanf(line,"Argument of Perigee(rad):   %f\n");
    
    line=fgetl(fid);
    dataset_out(i).mean_anomaly=sscanf(line,"Mean Anom(rad):             %f\n");
    
    line=fgetl(fid);
    dataset_out(i).af0=sscanf(line,"Af0(s):                    %f\n");
    
    line=fgetl(fid);
    dataset_out(i).af1=sscanf(line,"Af1(s/s):                  %f\n");
    
    line=fgetl(fid);
    dataset_out(i).week=sscanf(line,"week:                        %d\n");
    
    
    %Skip the last line and the first line of the next block
    line=fgetl(fid);
    line=fgetl(fid);
    
end