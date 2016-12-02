%This Matlab code was written by Elisabeth Gallant for use in my MS thesis
%please cite: 
%
%Gallant, E. (2016) MS Thesis, Advisor Charles Connor: Lava Flow Hazard Assessment for the Idaho National Laboratory, Idaho Falls, and Pocatello, Idaho, USA. University of South Florida Scholar Commons.

clear all;

% initialize colormap
cmap = jet; close;

% spatial distance template for event semimajor and minor axes
major = 5000;
minor = 2500;
NA = 361; %samples every degree (NA = number of angles)
orientation = 150*pi/180;

% set a temporal heirarchical clustering cutoff metric cutoff<(data-mean)/std
cutoff = 1.1;

% read in the data
Vage = xlsread('vents.xls','Age');
Vno = xlsread('vents.xls','No Age');

filename = ['temporal cluster major=' num2str(major) ' minor=' num2str(minor) ' orientation=' num2str(orientation) '.txt'];
fid = fopen(filename,'w');
[angles, template] = get_ellipse_template(major,minor,orientation,NA);

% extract the east and north locations
ENAge = Vage(:,4:5);
ENNo = Vno(:,4:5);

% extract the time for each vent
T = Vage(:,9);

% construct a heirarchical tree describing the time data
Z = linkage(T);

% perform temporal clustering
Tc = cluster(Z,'cutoff',cutoff);

NC = max(Tc);
mu = zeros(NC,1);
ENTime = cell(NC,1);
% calculate the mean of each cluster and store the members of each cluster
for k=1:NC
    curr = Tc==k;
    mu(k) = mean(T(curr));
    ENTime{k} = [ENAge(curr,:) T(curr)];
end

% rearrange the clusters so they are in chronological order
[mu, inx] = sort(mu);
Dir = ones(NC,1)*NaN;
Cluster = cell(NC,1);
Ndata = zeros(NC,1);

% get limits for consistant plotting (1000 chosen for astetic purposes)
dim = [min(ENAge(:,1))-1000, max(ENAge(:,1))+1000, min(ENAge(:,2))-1000, max(ENAge(:,2))+1000];
total_class = 0;
for k=1:max(Tc)
    Cluster{k} = ENTime{inx(k)};
    Ndata(k) = size(Cluster{k},1);
    % calculate the orientation if the cluster has more than one member
    if(Ndata(k)>1)       
        % perfrom spatial clustering
        not_in = true;
        X = Cluster{k}(:,1:2);
        Zc = linkage(X);
        nspace = 1;
        
        % plot dendrogram for each age group
        % plot the cluster member locations
        close all;
        f = figure;
        dendrogram(Zc,0);
        xlabel('Class')
        ylabel('Age Difference (years)')
        title(['Dendrogram for Temporal Cluster ' num2str(k) ' Time = ' num2str(mu(k))]);
        
        % save the image to a jpg
        saveas(f,['Dendrogram ' num2str(k) ' Time = ' num2str(mu(k)) '.jpg']); 
       
        while(not_in)
            % assume nspace clusters
            Ts = cluster(Zc,'maxclust',nspace);

            % loop through each spatial cluster
            ms = zeros(nspace,2);
            violated = false;
            for j=1:nspace
                %extract members of the spatial cluster
                member = X(Ts==j,:);
                
                % number of members
                Nm = size(member,1);
                                
                % find the mean of the cluster
                if(Nm>1)
                    ms(j,:) = mean(member);
                    
                    % compute vector from cluster mean to member
                    vec = member-ms(j*ones(Nm,1),:);
                    
                    % calculate angle between north and cluster-to-member
                    % vectors
                    A = atan2(vec(:,1),vec(:,2))*180/pi;
                    A(A<0) = A(A<0)+360;
                    
                    % calculate distance to mean
                    dist = sqrt(vec(:,1).^2+vec(:,2).^2);
                    
                    % obtain maximum template distance
                    Tdist = interp1(angles,template,A);
                    
                    % if any member is out of the template bounds stop
                    % checking
                    if(any(dist>Tdist))
                        violated = true;
                        break;
                    end
                else
                    ms(j,:) = member;
                    continue;
                end
            end
            
            % if any members were out of their template bounds redo
            % clustering assuming another cluster
            if(violated)
                nspace = nspace+1;
            else
                not_in = false;
            end           
        end
        
        % append cluster info to include spatial cluster info
        Cluster{k}(:,4) = Ts+total_class;
        total_class = total_class+max(Ts);
    else
        total_class = total_class+1;
        Cluster{k}(4) = total_class;
        Ts = 1;
    end
    
    % add cluster values to file (easting, northing, time, 
    % spatial cluster member ID)
    fprintf(fid,'%12.4f %12.4f %12.8f %4i\n',Cluster{k}');
   
    % get color indices
    cinx = round(linspace(1,64,nspace));
    
    % plot the cluster member locations
    close all;
    h = figure;
    hold on;
    for j=1:nspace
        plot(Cluster{k}(Ts==j,1),Cluster{k}(Ts==j,2),'.','markeredgecolor',cmap(cinx(j),:));
    end
    xlabel('Easting (m)')
    ylabel('Northing (m)')
    axis(dim);
    axis equal;
    title(['Vent Data Temporal Cluster ' num2str(k) ' Time = ' num2str(mu(k))]);
    
    % save the image to a jpg
    saveas(h,['Vents Temporal Cluster ' num2str(k) ' Time = ' num2str(mu(k)) '.jpg']); 
end

fid = fclose(fid);
