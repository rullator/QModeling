function results = QModeling_SRTM2(Ct,Cr,times,k2_p,k2a_min,k2a_max,basisF,resampling,varargin)
% It is used as a low-level function, called by QModeling_preprocessTACs or
% QModeling_parametric_images
%---------------------------------------------------------------------------------
% QModeling_SRTM2 is part of QModeling.
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


    num_files=length(Ct);

    theta_3 = exp(log(k2a_min):(log(k2a_max)-log(k2a_min))/(basisF-1):log(k2a_max)); %rango logaritmico de valores para el theta_3        
    theta_3 = theta_3/60; 

    tseg = zeros(1,num_files);

    for i=1:num_files 
        tseg(i)=(times(i,1)+times(i,2))/2;
    end

    RSS_menor = Inf; 
    i_min = 1;
    W=eye(num_files); %weight diagonal matrix. Is the identity matrix (in fact, the weigths are not consider but we follow Gunn algorithm
    w = ones(1,num_files); %vector of weights initializing by one
    B=zeros(length(theta_3),num_files);

    tgrid = 0:resampling:times(num_files,2);
    Cr_grid = interp1([0 tseg],[0;Cr],tgrid,'linear');
    Cr_grid(isnan(Cr_grid)) = Cr_grid(find(isnan(Cr_grid),1,'first')-1); %all nonvalues of the vector Cr_grid are changed by the real number which appears in the vector before the first NaN

    if (k2_p == -1) %is selected so we must calculate and estimate k2'
        M=zeros(2*length(theta_3),num_files); %for each theta_3 value, fill in two rows

        for i=1:1:length(theta_3) %for each basis function

            conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid));    
            B(i,:) = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %basis function i calculated for each frame

            A = [Cr B(i,:)'];      
            [Q,R]=qr(W*A,0); %QR decomposition reduced
            M(2*i-1:2*i,:)=R\Q'; %(2 x num_frames)

            theta=M(2*i-1:2*i,:)*W*Ct; %(2 x 1)

            Ct_i=theta(1)*Cr + theta(2)*B(i,:)'; %Stimated concentration for the basis function B_i in each frame (num_frames x 1)

            L=(Ct-Ct_i).^2; %Matrix to the square error (num_frames x 1)

            RSS_aux = w*L;        

            if (RSS_aux < RSS_menor)  %Look for the basis function's index which produce the minimum error
                    RSS_menor = RSS_aux;                      
                    i_min = i;
            end  
        end

        clear A Q R RSS_menor RSS_aux;

        theta=(M(2*i_min-1:2*i_min,:)*W*Ct);    
        theta_1=theta(1);   
        theta_2=theta(2);    

        %Obtain the parameters
        theta_3_fin=theta_3(i_min);

        R1_fin=theta_1;    
        k2a_fin=theta_3_fin;    
        k2_fin=theta_1*k2a_fin+theta_2;    
        BPnd_fin=(k2_fin/k2a_fin)-1;    
        k2_P = (k2_fin/R1_fin);  

        Ct_estimate=theta(1)*Cr+theta(2)*B(i_min,:)' ; %stimated concentration for each frame
        
        % Goodness of fit
        NMSE =1-((Ct-Ct_estimate).^2)./((Ct-mean(Ct)).^2);
        CCoef=corrcoef(Ct,Ct_estimate);

        % obtain the basis function using the k2' calculated
        for i=1:1:length(theta_3) %for each basis function
            conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid));    
            B(i,:) = Cr' + (k2_P-theta_3(i))*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %basis function i calculated for each frame
        end

    else %k2' is not selected, so we consider the given value
          k2_P=k2_p/60;    

       for i=1:1:length(theta_3) %for each basis function

            conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid));    
            B(i,:) = Cr' + (k2_P-theta_3(i))*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %basis function i calculated for each frame

            RI = (B(i,:)*Ct)/sum(B(i,:).^2);

            L=(Ct-RI*B(i,:)').^2; %Matrix to the square error (num_frames x 1)

            RSS_aux = w*L;
            if RSS_aux < RSS_menor %Look for the basis function's index which produce the minimum error
                RSS_menor = RSS_aux;
                i_min = i;
            end 
        end

            %Obtain and save the stimated parameters
            theta_3_fin=theta_3(i_min); 

            R1_fin=(B(i_min,:)*Ct)/sum(B(i_min,:).^2);
            k2a_fin=theta_3_fin;
            BPnd_fin=R1_fin*(k2_P/k2a_fin)-1;

            Ct_estimate=R1_fin*B(i_min,:); %stimated concentration for each frame
       
            % Goodness of fit
            % MSE=mean((Ct-Ct_estimate').^2);
            NMSE =1-((Ct-Ct_estimate').^2)./((Ct-mean(Ct)).^2);
            CCoef=corrcoef(Ct,Ct_estimate');
    end  

    if nargin==8 % preprocess
%             results = [R1_fin k2a_fin*60 k2_P*60 BPnd_fin MSE CCoef(1,2)];
              results = [R1_fin k2a_fin*60 k2_P*60 BPnd_fin mean(NMSE) CCoef(1,2)];
    elseif nargin==9 % parametric images
        Cpet=varargin{1};
        num_vox=size(Cpet,2);
        k2_p=k2_P*60;
        
        clear W w i_min theta Ct_i RI

        RSS_menor=Inf(1,num_vox); %Initialize the minimum square error for the basis function i(num_vox x 1) 
        i_min=ones(num_vox,1);
        RI_temp=Inf(1,num_vox);
        w=ones(1,num_files);

        for i=1:1:length(theta_3) %for each basis function
            RI = (B(i,:)*Cpet)/sum(B(i,:).^2);    
            L=(Cpet-B(i,:)'*RI).^2; %Matrix to the square error (num_frames x 1)
            RSS_aux = w*L;

            indices_men=find(RSS_aux<RSS_menor); 
            %Take the minimum values of the mean square error for each voxel and the basis function's number which generate it
            RSS_menor(indices_men)=RSS_aux(indices_men);
            i_min(indices_men)=i;

            RI_temp(indices_men)=RI(indices_men); 
        end

        theta_3_fin=theta_3(i_min); %for each voxel, load the value of theta_3 which produces the minimum square error 

        R1=RI_temp;
        k2a=theta_3_fin*60;
        %k2a=theta_3_fin;
        k2=k2_p*R1;
        BPnd=((R1./k2a)*k2_p) -1;

        %Save the images
        images(1,:)=BPnd;
        images(2,:)=k2;
        images(3,:)=k2a;
        images(4,:)=R1;

        results=images;
    end
end
