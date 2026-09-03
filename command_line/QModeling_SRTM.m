function results = QModeling_SRTM(Ct,Cr,times,k2a_min,k2a_max,num_BF,resampling,varargin)
% It is used as a low-level function, called by QModeling_preprocessTACs or
% QModeling_parametric_images
%---------------------------------------------------------------------------------
% QModeling_SRTM is part of QModeling.
%
% QModeling is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QModeling is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QModeling.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------------------------
% Copyright (C) 2016 Francisco Javier Lopez Gonzalez, Karl-Khader
% Thurnhofer Hemsi


    %Obtain info
    num_files = length(Ct);
    t = ((times(1:num_files,1)+times(1:num_files,2))/2)/60;

    %Logarithmic scale
    theta_3 = exp(log(k2a_min):(log(k2a_max)-log(k2a_min))/(num_BF-1):log(k2a_max)); 
    theta_3 = theta_3/60;

    %Create auxiliary matrix
    w = ones(1,num_files); 
    B=zeros(length(theta_3),num_files);
    M=zeros(2*length(theta_3),num_files); 
    W=eye(num_files); 

    tseg = t'*60;
    RSS_menor = Inf; 
    i_min = 1;

    % Interpolate
    tgrid = 0:resampling:times(end,2);
    Cr_grid = interp1([0 tseg],[0;Cr],tgrid,'linear');
    Cr_grid(isnan(Cr_grid)) = Cr_grid(find(isnan(Cr_grid),1,'first')-1);

    %Calculations for the algorithm
    for i=1:1:length(theta_3) 

        conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid)); %Convolution
        B(i,:) = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate

        A = [Cr B(i,:)'];      
        [Q,R]=qr(W*A,0); % QR factoritation
        M(2*i-1:2*i,:)=R\Q'; 

        theta=M(2*i-1:2*i,:)*W*Ct; 

        Ct_i=theta(1)*Cr+theta(2)*B(i,:)'; 

        L=(Ct-Ct_i).^2; 

        RSS_aux = w*L;
        if RSS_aux < RSS_menor
            RSS_menor = RSS_aux;
            i_min = i;
        end  
    end


    clear A Q R RSS_menor RSS_aux;

    theta=(M(2*i_min-1:2*i_min,:)*W*Ct);
    theta_1=theta(1);    
    theta_2=theta(2);    

    Ct_estimate=theta_1*Cr+theta_2*B(i_min,:)'; 
    

    % Get the parameters
    theta_3_fin=theta_3(i_min);
    R1_fin=theta_1;
    k2_fin=60*(theta_1*theta_3_fin+theta_2);
    k2a_fin=theta_3_fin*60;
    BPnd_fin=(k2_fin/k2a_fin)-1;

    % Goodness of fit
%     MSE=mean((Ct-Ct_estimate).^2);
    NMSE =1-((Ct-Ct_estimate).^2)./((Ct-mean(Ct)).^2);
    CCoef=corrcoef(Ct,Ct_estimate);

    if nargin==7 % preprocess
         results = [R1_fin k2a_fin k2_fin BPnd_fin mean(NMSE) CCoef(1,2)];
    elseif nargin==8 % parametric images
        Cpet=varargin{1};
        num_vox=size(Cpet,2);
        
        clear W w i_min theta Ct_i

        W=eye(num_files);
        w = ones(1,num_files);

        RSS_menor=Inf(1,num_vox); 
        i_min=ones(num_vox,1);
        theta_1_fin=Inf(num_vox,1);
        theta_2_fin=Inf(num_vox,1);


        for i=1:1:length(theta_3)

            theta=M(2*i-1:2*i,:)*W*Cpet;

            A = [Cr B(i,:)'];
            Ct_i = A*theta;

            RSS_aux=(w*((Cpet-Ct_i).^2));
            clear Ct_i;

            indices_men=find(RSS_aux<RSS_menor); 

            RSS_menor(indices_men)=RSS_aux(indices_men);
            i_min(indices_men)=i;

            theta_1_fin(indices_men)=theta(1,indices_men);
            theta_2_fin(indices_men)=theta(2,indices_men);

        end
        
        theta_3_fin=theta_3(i_min); 

        % Get parametric images
        R1=theta_1_fin;
        k2=60*(theta_1_fin'.*theta_3_fin+theta_2_fin');
        k2a=theta_3_fin*60;
        BPnd=(k2-k2a)./k2a;

        images(1,:)=BPnd;    
        images(2,:)=k2;
        images(3,:)=k2a;
        images(4,:)=R1;
        
        results=images;

    end
end
