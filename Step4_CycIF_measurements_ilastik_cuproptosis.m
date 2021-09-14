function Step4_CycIF_measurements_ilastik_cuproptosis(filename, options)  
    rng(1);

    for k = 1:length(filename.tissues)
        fold = filename.tissues{k};
        filename.folders.resultfile = [filename.folders.results fold '_Results_' options.date '.mat'];

        addpath(filename.folders.main)

        %%%% initialize the variables we are going to measure

            Results = [];
    %         Results.FieldCoord = [];
    %         Results.FieldFlag = [];
            Results.Area = [];
            Results.Solidity = [];
            Results.MedianNucSign = [];
            Results.MedianCytSign = [];
            Results.MeanNucSign = [];
            Results.MeanCytSign = [];
            Results.TotalNucSign = [];
            Results.TotalCytSign = [];
            Results.MedianNucRaw = [];
            Results.MedianCytRaw = [];
            Results.MeanNucRaw = [];
            Results.MeanCytRaw = [];
            Results.TotalNucRaw = [];
            Results.TotalCytRaw = [];

            Results.RepIndex = [];
            Results.Foci.Area = [];
            Results.Foci.Centroids = [];


        reps = 'abc';
        %%%% loop around the split fields
        for r = reps
            for f = 1:56
    %         for f = options.randfields
                tic
                %%%% load the FullStack and/to check it exists ----------------
                FieldFile = [fold r '_Field_' num2str(f)];
                FileTif = [filename.folders.main filename.folders.composites FieldFile  filename.suffix];

                if exist(FileTif,'file')~=2
                    disp(FileTif)
                    disp('Error: FullStack not found')
                    continue
                end
                disp(FieldFile)

                %%%% load nuclear, cytoplasmic, and foci masks
                segfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg FieldFile '_Seg.tif'];
                cyt_segfile = [filename.folders.main filename.folders.output filename.folders.ilastikseg FieldFile '_CytSeg.tif'];
                foci_segfile = [filename.folders.main filename.folders.output 'FociSegmentation_20210721\' FieldFile '_FociSeg.tif'];
                assignedfoci_segfile = [filename.folders.main filename.folders.output 'FociSegmentation_20210721\' FieldFile '_AssignedFociSeg.tif'];

                if exist(segfile,'file') ~= 2
                    disp('Nuclear segmentation file not found')
                    continue
                end

                if exist(cyt_segfile,'file') ~= 2
                    disp('Cytoplasmic segmentation file not found')
                    continue
                end

                NucMask = uint16(imread(segfile));
                lb_NucMask = NucMask;

                % load cytoplasmic mask
                CytMask = uint16(imread(cyt_segfile));
                lb_CytMask = CytMask;

                fociflag = 0;
                if (exist(foci_segfile,'file') == 2 && exist(assignedfoci_segfile,'file') == 2)
                   lb_foci = uint16(imread(foci_segfile)); 
                   assigned_foci = uint16(imread(assignedfoci_segfile));
                   fociflag = 1;
                end

                %%%% measure cell properties and assign them to results -------

                stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
                stats_cyt = regionprops(lb_CytMask,'Area');

                if(fociflag)
                   stats_foci = regionprops(lb_foci,'Area','Centroid'); 
                end

                Area = {stats_nuc.Area};
                Solidity = {stats_nuc.Solidity};
                Centroid = {stats_nuc.Centroid};
                CentroidVect = cell2mat(Centroid);
                CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
                CytArea = {stats_cyt.Area};

                totcells_field = length(Area);
                currentcells = length(Results.Area);
                startcell = currentcells + 1;
                endcell = currentcells + totcells_field;

                %foci stats
                if(fociflag)
                   stats_foci = regionprops(lb_foci,'Area','Centroid'); 

                   FociArea = {stats_foci.Area};
                   FociCentroid = {stats_foci.Centroid};
                   FociCentroidVect = cell2mat(FociCentroid);
                   FociCentroidMat = reshape(FociCentroidVect,2,length(FociCentroidVect)/2);

                   totfoci = length(unique(lb_foci))-1;

                   fociidx = zeros(length(totfoci),1);

                   counter = 1;
                   for i = unique(lb_foci)'
                       if i == 0 
                           continue
                       end
                      fociidx(counter) = assigned_foci(ceil(FociCentroidMat(2,i)),ceil(FociCentroidMat(1,i))) + currentcells; 
                      counter = counter + 1;
                   end

                   currentfoci = length(Results.Foci.Area);
                   startfocus = currentfoci + 1;
                   endfocus = currentfoci + totfoci;

                   Results.Foci.Area(startfocus:endfocus,1) = cell2mat(FociArea)';
                   Results.Foci.Centroids(startfocus:endfocus,1) = ceil(FociCentroidMat(1,:));
                   Results.Foci.Centroids(startfocus:endfocus,2) = ceil(FociCentroidMat(2,:));


                end


                if totcells_field < 1
                    disp(['No cells in Field ' FieldFile])
                    continue
                end

                Results.Area(startcell:endcell,1) = cell2mat(Area)';
                Results.Solidity(startcell:endcell,1) = cell2mat(Solidity);
                Results.CentroidX(startcell:endcell,1) = CentroidMat(1,:); %+Coordinates.Field{i1,i2}(4);
                Results.CentroidY(startcell:endcell,1) = CentroidMat(2,:); %+Coordinates.Field{i1,i2}(2);
                Results.CytArea(startcell:endcell,1) = [cell2mat(CytArea)'; zeros(totcells_field-length(cell2mat(CytArea)),1)];
    %             Results.FieldCoord(startcell:endcell,1) = zeros(totcells_field,1)+i1;
    %             Results.FieldCoord(startcell:endcell,2) = zeros(totcells_field,1)+i2;


                Results.RepIndex(startcell:endcell,1) = find(reps == r);

                %%%% cycle around the images and extract measurements ---------
                stackinfo = imfinfo(FileTif);
                maxround = length(stackinfo);
                for i3 = 1:maxround
                    % load image
                    FluorImage = imread(FileTif,'Index',i3);
                    FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
                    % extract pixel values from fluo image using the segmented images    
                    Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
                    Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');

                    Nuclei_stats_raw = regionprops(lb_NucMask,FluorImage,'PixelValues');
                    Cytopl_stats_raw = regionprops(lb_CytMask,FluorImage,'PixelValues');

                    %pad arrays
                    temp = nan(endcell-startcell+1,1);
                    Results.MedianNucSign(startcell:endcell,i3) = temp;
                    Results.MeanNucSign(startcell:endcell,i3) = temp;
                    Results.TotalNucSign(startcell:endcell,i3) = temp;

                    %add signal values    
                    nuc_end = startcell + length({Nuclei_stats.PixelValues})-1;
                    cyt_end = startcell + length({Cytopl_stats.PixelValues})-1;

                    Results.MedianNucSign(startcell:nuc_end,i3) = cellfun(@median,{Nuclei_stats.PixelValues});
                    Results.MeanNucSign(startcell:nuc_end,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
                    Results.TotalNucSign(startcell:nuc_end,i3) = cellfun(@sum,{Nuclei_stats.PixelValues});

                    Results.MedianNucRaw(startcell:nuc_end,i3) = cellfun(@median,{Nuclei_stats_raw.PixelValues});
                    Results.MeanNucRaw(startcell:nuc_end,i3) = cellfun(@mean,{Nuclei_stats_raw.PixelValues});
                    Results.TotalNucRaw(startcell:nuc_end,i3) = cellfun(@sum,{Nuclei_stats_raw.PixelValues});


                    Results.MedianCytSign(startcell:cyt_end,i3) = cellfun(@median,{Cytopl_stats.PixelValues});
                    Results.MeanCytSign(startcell:cyt_end,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
                    Results.TotalCytSign(startcell:cyt_end,i3) = cellfun(@sum,{Cytopl_stats.PixelValues});

                    Results.MedianCytRaw(startcell:cyt_end,i3) = cellfun(@median,{Cytopl_stats_raw.PixelValues});
                    Results.MeanCytRaw(startcell:cyt_end,i3) = cellfun(@mean,{Cytopl_stats_raw.PixelValues});
                    Results.TotalCytRaw(startcell:cyt_end,i3) = cellfun(@sum,{Cytopl_stats_raw.PixelValues});                

                end

                %foci measurements
                if(fociflag)
                    FluorImage = imread(FileTif,'Index',2);
                    FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));

                    Foci_stats = regionprops(lb_foci,FluorImage_BK,'PixelValues');
                    Foci_stats_raw = regionprops(lb_foci,FluorImage,'PixelValues');

                    focimean = cellfun(@median,{Foci_stats.PixelValues});
                    focimedian = cellfun(@mean,{Foci_stats.PixelValues});
                    totsign = cellfun(@sum,{Foci_stats.PixelValues});                

                    Results.Foci.MedianSign(startfocus:endfocus,1) = focimean(~isnan(focimean));
                    Results.Foci.MeanSign(startfocus:endfocus,1) = focimedian(~isnan(focimedian));
                    Results.Foci.TotalSign(startfocus:endfocus,1) = totsign(~isnan(totsign));

                    focimeanraw = cellfun(@median,{Foci_stats_raw.PixelValues});
                    focimedianraw = cellfun(@mean,{Foci_stats_raw.PixelValues});
                    totsignraw = cellfun(@sum,{Foci_stats_raw.PixelValues});                

                    Results.Foci.MedianRaw(startfocus:endfocus,1) = focimeanraw(~isnan(focimeanraw));
                    Results.Foci.MeanRaw(startfocus:endfocus,1) = focimedianraw(~isnan(focimedianraw));
                    Results.Foci.TotalRaw(startfocus:endfocus,1) = totsignraw(~isnan(totsignraw));               

                    Results.Foci.CellID(startfocus:endfocus,1) = fociidx;
                    Results.Foci.WellID(startfocus:endfocus,1) = k;
                    Results.Foci.RepID(startfocus:endfocus,1) = r;
                    Results.Foci.FieldID(startfocus:endfocus,1) = f;
                end

    %             Results.FieldFlag(i1,i2)=1; 

                toc
                disp(['Finished ' FieldFile])
            end

        end

         save([filename.folders.main filename.folders.output filename.folders.resultfile],'Results', '-v7.3') %Saving matrix   

    end
end