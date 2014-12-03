classdef CellNetwork<handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        BW=20e6;%Hz
        N=2048;%Resource Elements
        Wk=15000;%Hz
        PTx=100;%Watts
        distance=2000;%diameter meters
        alpha=2.2;
        P_d=(10/1000)^2.5;%Path Loss at 10 metes
        Hk=1;%AWGN Channel
        No=1e-8;%Watts/Hz
        Nvar=[];
        Rb=30000000;%Bis/Sec
        SNR=[];
        SNIR=[];
        G_snr=[];
        G_snir=[];
        target_matrix_x=[];
        target_matrix_y=[];
        users=[];
    end
    
    methods
        function obj=CellNetwork(varargin)
            for k=1:nargin
                switch k
                    case 1
                        obj.BW=varargin{1};
                    case 2
                        obj.N=varargin{2};
                    case 3
                        obj.Wk=varargin{3};
                    case 4
                        obj.PTx=varargin{4};
                    case 5
                        obj.distance=varargin{5};
                    case 6
                        obj.alpha=varargin{6};
                    case 7
                        obj.P_d=varargin{7};
                    case 8
                        obj.Hk=varargin{8};
                    case 9
                        obj.No=varargin{9};
                    case 10
                        obj.Nvar=varargin{10};
                    case 11
                        obj.Rb=varargin{11};
                end
            end
            %function hexagon(cote,x0,y0)
            x0(1)=0.5*(obj.distance);
            y0(1)=0.5*(obj.distance);
            cote =0.5*(obj.distance);
            
            height=(sqrt(3)/2)*cote;
            for i=1:3
                x0(i+1)=x0(1);
                y0(i+1)=y0(1)+(i*2*height);
            end
            
            x0(5)=x0(1)+2*cote-(cote/2);
            y0(5)=y0(1)-height;
            
            for i=1:4
                x0(i+5)=x0(5);
                y0(i+5)=y0(5)+(i*2*height);
            end
            
            x0(10:13)=x0(1)+4*cote-(2*cote/2);
            x0(14:18)=x0(1)+6*cote-(3*cote/2);
            x0(19:22)=x0(1)+8*cote-(4*cote/2);
            x0(23:27)=x0(1)+10*cote-(5*cote/2);
            
            y0(10:13)=y0(1:4);
            y0(14:18)=y0(5:9);
            
            y0(19:22)=y0(1:4);
            y0(23:27)=y0(5:9);
            
            x0h=x0(:);
            y0h=y0(:);
            
            k1=1;
            
            distance_arr=0:10:cote;
            distance_arr(1)=1;
            
            figure()
            hold on
            for k=[1,2,5,6,7,10,11]
                xp=x0(k);
                yp=y0(k);
                x(k,:)=cote*[-1 -0.5 0.5 1 0.5 -0.5 -1]+xp;
                y(k,:)=cote*sqrt(3)*[0 -0.5 -0.5 0 0.5 0.5 0]+yp;
                x_current=x(k,:);
                y_current=y(k,:);
                plot(x_current,y_current,'b','Linewidth',2);%grid;
                xlabel('X Label meters')
                ylabel('Y Label meters')
            end
            
            Nai=(1.38e-23)*(((10^(4/10))*(10^(8/10)))-1)*(290)*obj.BW;
            Nin=(1.38e-23)*(1250)*obj.BW;
            obj.Nvar=((10^(11/10))/(10^(4/10)))*(Nai+Nin);
            temp_PTx=(obj.PTx*((10^(11/10))/(10^(4/10))))/obj.Rb;
            obj.PTx=temp_PTx;
            temp_No=obj.Nvar/obj.BW;
            obj.No=temp_No;
            for t=0:359
                targetx=x0(6);
                targety=y0(6);
                %target_matrix(t+1,:)=((t*pi)/180)*(distance_arr+targetx);%cosd(t)
                target_matrix_x(t+1,:)=(distance_arr*cosd(t))+targetx;
                target_matrix_y(t+1,:)=(distance_arr*sind(t))+targety;
                distance_node1(t+1,:)=sqrt(((target_matrix_y(t+1,:)-y0(1)).^2)+((target_matrix_x(t+1,:)-x0(1)).^2));
                distance_node2(t+1,:)=sqrt(((target_matrix_y(t+1,:)-y0(2)).^2)+((target_matrix_x(t+1,:)-x0(2)).^2));
                distance_node5(t+1,:)=sqrt(((target_matrix_y(t+1,:)-y0(5)).^2)+((target_matrix_x(t+1,:)-x0(5)).^2));
                distance_node7(t+1,:)=sqrt(((target_matrix_y(t+1,:)-y0(7)).^2)+((target_matrix_x(t+1,:)-x0(7)).^2));
                distance_node10(t+1,:)=sqrt(((target_matrix_y(t+1,:)-y0(10)).^2)+((target_matrix_x(t+1,:)-x0(10)).^2));
                distance_node11(t+1,:)=sqrt(((target_matrix_y(t+1,:)-y0(11)).^2)+((target_matrix_x(t+1,:)-x0(11)).^2));
                interference_node1(t+1,:)=(obj.PTx*(distance_node1(t+1,:).^(-obj.alpha)));
                interference_node2(t+1,:)=(obj.PTx*(distance_node2(t+1,:).^(-obj.alpha)));
                interference_node5(t+1,:)=(obj.PTx*(distance_node5(t+1,:).^(-obj.alpha)));
                interference_node7(t+1,:)=(obj.PTx*(distance_node7(t+1,:).^(-obj.alpha)));
                interference_node10(t+1,:)=(obj.PTx*(distance_node10(t+1,:).^(-obj.alpha)));
                interference_node11(t+1,:)=(obj.PTx*(distance_node11(t+1,:).^(-obj.alpha)));
                interference_sum(t+1,:)=interference_node1(t+1,:)+interference_node2(t+1,:)+interference_node5(t+1,:)+...
                    interference_node7(t+1,:)+interference_node10(t+1,:)+interference_node11(t+1,:);
                
                SNIR(t+1,:)=((obj.PTx*(distance_arr(1,:).^(-obj.alpha)))./((obj.PTx*(distance_node1(t+1,:).^(-obj.alpha)))+(obj.PTx*(distance_node2(t+1,:).^(-obj.alpha)))+...
                    (obj.PTx*(distance_node5(t+1,:).^(-obj.alpha)))+(obj.PTx*(distance_node7(t+1,:).^(-obj.alpha)))+...
                    (obj.PTx*(distance_node10(t+1,:).^(-obj.alpha)))+(obj.PTx*(distance_node11(t+1,:).^(-obj.alpha)))+obj.Nvar));
                
                G_snir(t+1,:)=10*log10(SNIR(t+1,:));
                
                SNR(t+1,:)=((obj.PTx*(distance_arr(1,:).^(-obj.alpha)))./(obj.Nvar));
                
                G_snr(t+1,:)=10*log10(SNR(t+1,:));
                
                users(t+1,:) = roundn(randi(1000,3,1),1);
                users_mapped_x(t+1,:)=(users(t+1,:)*cosd(t))+targetx;
                users_mapped_y(t+1,:)=(users(t+1,:)*sind(t))+targety;
                plot(users_mapped_x(t+1,:),users_mapped_y(t+1,:),'g.','markers',4);
            end
            obj.SNIR=SNIR;
            obj.SNR=SNR;
            obj.G_snir=G_snir;
            obj.G_snr=G_snr;
            obj.target_matrix_x=target_matrix_x;
            obj.target_matrix_y=target_matrix_y;
            obj.users=users;
            
        end
    end
    
end

