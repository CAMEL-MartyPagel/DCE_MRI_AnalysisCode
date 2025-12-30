function [DCE_images,TR,flip_angle] = load_DCE_image(DIR_DCE)

    % Obtain the number of echos and number of slices.
    cd(DIR_DCE);
    
    if exist([DIR_DCE,'/method'],'file')
        methodfile = [DIR_DCE,'method'];
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
            flag_RARE=isempty(strfind(sequence{1,1},'FLASH'));
            if  flag_RARE==1
                error('Error: The sequence selected is not a FLASH.');
            end
        end 
        if strncmp(current_line, '##$PVM_NEchoImages',18) 
            anatomical_number_of_echos = strread(current_line,'##$PVM_NEchoImages=%f');
        end
        if strncmp(current_line, '##$PVM_RepetitionTime',21) 
            TR = strread(current_line,'##$PVM_RepetitionTime=%f');
        end
        if strncmp(current_line, '##$PVM_NRepetitions',19)
            number_of_frames = strread(current_line,'##$PVM_NRepetitions=%f');
        end
        if strncmp(current_line, '##$PVM_SPackArrNSlices',22) 
            current_line=fgetl(method_file);
            anatomical_slices_packets=strread(current_line);
            anatomical_number_of_slices=anatomical_slices_packets(1);
        end
    end
    fclose(method_file);
    
    DIR_pdata=strcat(DIR_DCE,'/pdata/1/');
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
        if strncmp(current_line, '##$VisuAcqFlipAngle',19)
            flip_angle = strread(current_line,'##$VisuAcqFlipAngle=%f');
        end
    end
    fclose(visu_pars_file);
   
    % Set up a 3D matrix for the images
    DCE_images=zeros(anatomical_columns,anatomical_rows,anatomical_number_of_slices, number_of_frames);
    
    % Get the image data from the 2dseq file.
    cd(DIR_pdata);
    if exist([DIR_pdata,'/2dseq'],'file')
        imagefile = [DIR_pdata,'2dseq'];
    else
        error('Error! n\2dseq file does not exist, please check the path %s ',[DIR_pdata,'/2dseq'])
    end
    
    FileID = fopen(imagefile,'r','l');
    DCE_images_from_2dseq = fread(FileID,'int16');
    fclose(FileID);
    
    DCE_images=reshape(DCE_images_from_2dseq,anatomical_columns,anatomical_rows,anatomical_number_of_slices, number_of_frames);
    DCE_images=permute(DCE_images,[2 1 3 4]);
end