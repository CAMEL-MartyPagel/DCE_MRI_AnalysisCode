function [anatomical_images,anatomical_offsets,multiple_TRs]=load_anatomical_image(DIR_anatomical)

    % Obtain the number of echos and number of slices.
    cd(DIR_anatomical);
    
    if exist([DIR_anatomical,'/method'],'file')
        methodfile = [DIR_anatomical,'method'];
    else
        error('Error! \nmethod file does not exist, please check the path %s ',[DIR,'method'])
    end
    
    method_file=fopen(methodfile,'r');
    line_number=0;
    while 1
        line_number=line_number+1;
        current_line = fgetl(method_file);
        if ~ischar(current_line),  break,   end
        if strncmp(current_line, '##$Method',9)
            sequence = textscan(current_line, '##$Method=%s');
            flag_RARE=isempty(strfind(sequence{1,1},'RARE'));
            if  flag_RARE==1
                error('Error: The sequence selected is not a RARE.');
            end
        end 
        if strncmp(current_line, '##$PVM_NEchoImages',18) 
            anatomical_number_of_echos = strread(current_line,'##$PVM_NEchoImages=%f');
        end
        if strncmp(current_line, '##$MultiRepTime',15)
            current_line=fgetl(method_file);
            multiple_TRs=strread(current_line);
            %number_of_TRs = strread(current_line,'##$$MultiRepTime=( %f )');
        end
        if strncmp(current_line, '##$T1Exp',8) 
            number_of_T1 = strread(current_line,'##$T1Exp=%f');
        end
        if strncmp(current_line, '##$PVM_SPackArrNSlices',22) 
            current_line=fgetl(method_file);
            anatomical_slices_packets=strread(current_line);
            anatomical_number_of_slices=anatomical_slices_packets(1);
        end
    end
    fclose(method_file);
    
    DIR_pdata=strcat(DIR_anatomical,'/pdata/1/');
    cd(DIR_pdata)
    
    if exist([DIR_pdata,'/visu_pars'],'file')
        visu_parsfile = [DIR_pdata,'visu_pars'];
    else
        error('Error! n\visu_pars file does not exist, please check the path %s ',[DIR_pdata,'/visu_pars'])
    end
      
    visu_pars_file=fopen(visu_parsfile,'r');
    line_number=0;
    while 1
        line_number=line_number+1;
        current_line = fgetl(visu_pars_file);
        if ~ischar(current_line),   break,   end
        if strncmp(current_line, '##$VisuCoreSize',15) 
            current_line=fgetl(visu_pars_file);
            anatomical_image_size=strread(current_line);
            anatomical_rows=anatomical_image_size(2);
            anatomical_columns=anatomical_image_size(1);
        end
        if strncmp(current_line, '##$VisuCorePosition',19)
            anatomical_position_linestart = line_number+1;
            current_line=fgetl(visu_pars_file);
            anatomical_offsets=strread(current_line);
        end
    end
    fclose(visu_pars_file);
    
    % Set up a 3D matrix for the images
    anatomical_images=zeros(anatomical_columns,anatomical_rows,anatomical_number_of_slices, number_of_T1);
    
    % Get the image data from the 2dseq file.
    cd(DIR_pdata);
    if exist([DIR_pdata,'/2dseq'],'file')
        imagefile = [DIR_pdata,'2dseq'];
    else
        error('Error! n\2dseq file does not exist, please check the path %s ',[DIR_pdata,'/2dseq'])
    end
    
    FileID = fopen(imagefile,'r','l');
    anatomical_images_from_2dseq = fread(FileID,'int16');
    fclose(FileID);
    
    anatomical_images=reshape(anatomical_images_from_2dseq,anatomical_columns,anatomical_rows,anatomical_number_of_slices, number_of_T1);
    anatomical_images=permute(anatomical_images,[2 1 3 4]);
end
