function Step3_segmentfromilastik_cuproptosis(filename, options) 
% open the tiff file of the full stack images and segments them based on
% Ilastik Probabilities file 
% Inputs FullStacks and Ilastik Probabilities files 

tic
mkdir([filename.folders.main filename.folders.output filename.folders.ilastikseg])

% for k = 1:length(filename.tissues) 
for k = 1:8
    fold = filename.tissues{k}; 
    addpath(filename.folders.main)
    reps = 'abc';

    for r = reps
        for f = 1:56
%          for f = options.randfields
            core = [fold r '_Field_' num2str(f)];
            FileTif = [filename.folders.main filename.folders.output  filename.folders.ilastikstacks core filename.suffix];
            
            segfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg core '_Seg.tif'];
            cytsegfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg core '_CytSeg.tif'];
%             focisegfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg core '_FociSeg.tif'];
%             assignedfocifile = [filename.folders.main filename.folders.output filename.folders.ilastikseg core '_AssignedFociSeg.tif'];
            focisegfile = [filename.folders.main filename.folders.output 'FociSegmentation_20210721\' core '_FociSeg.tif'];
            assignedfocifile = [filename.folders.main filename.folders.output 'FociSegmentation_20210721\' core '_AssignedFociSeg.tif'];

            IlastikTif = [filename.folders.main filename.folders.output filename.folders.ilastikprob core  filename.ilastiksuffix];
            FullTif = [filename.folders.main filename.folders.composites core filename.suffix];
            
            if (exist(segfile,'file')==2 && exist(cytsegfile,'file')==2 && exist(focisegfile,'file')==2)
                disp(['Image ' core ' already segmented'])
                continue
            elseif (exist(segfile,'file')==2 && exist(cytsegfile,'file')==2 && k > 16)
                disp(['Image ' core ' already segmented, and no foci segmentation requested'])
                continue
            elseif exist(IlastikTif,'file')~=2
                disp(IlastikTif)
                disp(['No ilastik prob file for ' core])
                continue
            else
                
            end
            
            disp(FileTif)

            try 
                DAPI = imread(FileTif,'Index',1);
                %ominfo = imfinfo(FileTif);
                %DAPI_last = imread(FileTif,'Index',length(ominfo)-3);
            catch
                disp('Error: DAPI not found')
                continue
            end
            
            
            try
                IlastikRGB = imread(IlastikTif);
            catch
                disp('Error: Ilastik RGB not found')
                disp(IlastikTif)
                continue
            end
            
            if(exist(segfile,'file') == 2)
                Nuclei = imread(segfile);
                disp('Using already segmented nuclei')
                if(exist(cytsegfile,'file') == 2)
                    disp('Using already segmented cytoplasm')
                    Cytoplasm = imread(cytsegfile);
                end
            else                          
                NucProb = IlastikRGB(:,:,options.nuc);
                BkgProb = IlastikRGB(:,:,options.backgr);            
                CytProb = IlastikRGB(:,:,options.cyt);

                NucProb = imopen(NucProb,strel('octagon',3));
                BkgProb = imopen(BkgProb,strel('octagon',3));
                                        
                % pre-filter easy ones
                nuc_mask = zeros(size(NucProb));
                temp_mask = zeros(size(NucProb));
                else_mask = zeros(size(NucProb));

                index = BkgProb > options.max_prob*options.bkgdprob_min | CytProb > options.max_prob*options.cytprob_min; 
                else_mask(index) = 1;

                index = NucProb > options.max_prob*options.nucprob_min; 
                temp_mask(index) = 1;

                nuc_mask = temp_mask & ~else_mask;
                nuc_mask = bwareaopen(nuc_mask,options.cellsize*2,4);

                sImage = NucProb;
                ThresImage=nuc_mask;

                seeds = imclose(sImage,strel('disk',round(options.cellsize/4)));
                seeds = imgaussfilt(seeds,12,'FilterSize',round(options.cellsize)); %%changed 5 --> 12
                seeds = imregionalmax(seeds);
                seeds=seeds&ThresImage;
                seeds=imdilate(seeds,strel('disk',2));

                % threshold and watershed image
                L=watershed(imimposemin(-sImage,seeds));
                Nuclei= uint32(ThresImage&L);

                Nuclei = bwareaopen(Nuclei,options.cellsize*2,4);
                Nuclei = imfill(Nuclei,'holes');
                forcheck = uint16(Nuclei);
                Nuclei = bwlabel(uint32(Nuclei));
                Nuclei = imdilate(Nuclei,strel('disk',2));
            end
            
            if(unique(Nuclei)==0)
                disp('No Nuclei found in image')
                continue
            end
                       
            %cytoplasmic mask
            if(exist(cytsegfile,'file')~=2)
                %generate prob files if Nuclei were already segmented but
                %cytoplasm wasn't
                
                if(exist(segfile,'file') == 2)
                    NucProb = IlastikRGB(:,:,options.nuc);
                    BkgProb = IlastikRGB(:,:,options.backgr);            
                    CytProb = IlastikRGB(:,:,options.cyt);

                    NucProb = imopen(NucProb,strel('octagon',3));
                    BkgProb = imopen(BkgProb,strel('octagon',3));
                end
                
                cyt_mask = zeros(size(CytProb));
                temp_mask = zeros(size(NucProb));
                else_mask = zeros(size(NucProb));

                index = BkgProb > options.max_prob*options.bkgdprob_min | NucProb > options.max_prob*options.nucprob_min; 
                else_mask(index) = 1;

                index = CytProb > options.max_prob*options.cytprob_min; 
                temp_mask(index) = 1;

                cyt_mask = temp_mask & ~else_mask;
                cyt_mask = bwareaopen(cyt_mask,options.cellsize*2,4);

                [D,idx] = bwdist(Nuclei,'euclidean');

                idx2 = idx;

                for i = 1:size(idx2,1)
                   for j = 1:size(idx2,2)
                       temp = Nuclei(idx(i,j));
                       if(D(i,j) > 1000)
                           temp = 0;
                       end
                       idx2(i,j) = temp;

                   end
                end

                idx2(~cyt_mask) = 0;

                Cytoplasm = idx2;
            end
            
            if(exist(segfile,'file')~=2)
               check_img = cat(3,uint16(Nuclei),uint16(DAPI),uint16(forcheck));
               check_file = [filename.folders.main filename.folders.output filename.folders.ilastikseg core  '_checkSeg.tif'];
               imwrite(check_img,check_file) 
               
               % label the nuclei
               saveastiff(uint32(Nuclei),segfile)
            end
            
            if(exist(cytsegfile,'file')~=2)
                PHAL = imread(FileTif,'Index',2);
                cyt_img = cat(3,uint16(Nuclei),uint16(Cytoplasm),uint16(PHAL));
                check_cyt_file = [filename.folders.main filename.folders.output filename.folders.ilastikseg core  '_checkCytSeg.tif'];
                imwrite(cyt_img,check_cyt_file)
            
                saveastiff(uint32(Cytoplasm),cytsegfile)
            end
            
            %segment foci
            if (k >= 1 && k <= 16)
                disp('Starting foci segmentation')
                FOCI = uint16(imread(FullTif,'Index',2));
                FOCI_BK = FOCI-imopen(FOCI,strel('disk',round(options.focibkgsigma))); 

                spots = filterLoG(FOCI_BK, 20/options.focifiltfactor);
                
                if (k >= 1 && k <= 8)
                    threshold = options.dlatthresh;
                elseif(k >= 9 && k <= 16)
                    threshold = options.pdhathresh;
                end
                foci_seg = spots > threshold;

                foci_img = cat(3,uint16(foci_seg),uint16(FOCI_BK),uint16(Cytoplasm));
                check_foci_file = [filename.folders.main filename.folders.output 'FociSegmentation_20210721\' core  '_checkFociSeg.tif'];
                imwrite(foci_img,check_foci_file)

                foci_seg = bwlabel(foci_seg);
                saveastiff(uint32(foci_seg),focisegfile);
                
                %assign foci to nuclei
                
                %combine nuc and cyt masks
                cell_mask = Cytoplasm;
                for i = 1:length(unique(Nuclei))
                    cell_mask(Nuclei == i) = i;
                end      
                
                cell_mask = imfill(cell_mask,'holes');
                
                %dist transform of each focus from combined mask
                [D,idx] = bwdist(cell_mask,'euclidean');
                
                %get centroids of each focus to assign to a cell
                stats_foci = regionprops(foci_seg,'Centroid');

                Centroid = cell2mat({stats_foci.Centroid});
                CentroidMat = reshape(Centroid,2,length(Centroid)/2);

                CentroidX = CentroidMat(2,:)';
                CentroidY = CentroidMat(1,:)';
                
                assigned_foci = foci_seg;
                
                focis = unique(foci_seg);
                focis = focis(2:length(focis))';
                for i = focis
                    focusid = foci_seg(ceil(CentroidX(i)),ceil(CentroidY(i)));
                    index = idx(ceil(CentroidX(i)),ceil(CentroidY(i)));
                    cellid = cell_mask(index);
                    if(D(ceil(CentroidX(i)),ceil(CentroidY(i))) > 15)
                        cellid = 0;
                    end
                    assigned_foci(foci_seg == focusid) = cellid;
                end
                
                %check all masks combined
                check_mask = cell_mask;
            
                unique_nuclei = unique(Nuclei)';
                unique_nuclei = unique_nuclei(2:length(unique_nuclei));
                for i = unique_nuclei
                   check_mask(assigned_foci == i) = i; 
                end
                
                %save check image
                check_img = cat(3,uint16(foci_seg),uint16(Cytoplasm),uint16(check_mask));
                check_assign_file = [filename.folders.main filename.folders.output 'FociSegmentation_20210721\' core  '_checkAssignSeg.tif'];
                imwrite(check_img,check_assign_file)
                
                %save assignment
                saveastiff(uint32(assigned_foci),assignedfocifile);
            end
            
            if options.pausetime > 0
               disp(['Pausing for ' num2str(options.pausetime) ' seconds'])
               pause(options.pausetime)
            end
            
        end
    end
end 
