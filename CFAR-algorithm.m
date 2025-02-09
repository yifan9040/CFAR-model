clear
close all
clc
for noise=9
    clearvars -except noise
    disp(['exp:' num2str(noise)])
    %%%%%%%%%%%%%%%%%%%%load net%%%%%%%%%%%%%%%
    net1=load('LSTM_inter.mat');
    net2=load('LSTM_mean.mat');
    %%%%%%%%%%%%%%%%%%%%%%Experimental parameters%%%%%%%%%%%%%%%%%
    P_FA=1e-4;
    N_reference=18;
    SNR_dB=0:40;
    SNR=10.^(SNR_dB./10);
    M_num=1e4;
    ii=N_reference+1;
    T_CA=(P_FA.^(-1/(2*N_reference))-1);
    T_GOCA=(P_FA.^(-1/(N_reference))-1);
    T_SOCA=(P_FA.^(-1/(N_reference))-1);
    KVI=4.76;
    KMR=1.806;
    elapsedTime=0;
    %%%%%%%%%%%%%%%%%%%%%%initialization%%%%%%%%%%%%%%%%%
    ideal=zeros(1,length(SNR));
    detect_rate_CA=zeros(1,length(SNR));
    detect_rate_OS=zeros(1,length(SNR));
    detect_rate_VI=zeros(1,length(SNR));
    VI_P_CA_AB=zeros(1,length(SNR));
    VI_P_GO_AB=zeros(1,length(SNR));
    VI_P_CA_A=zeros(1,length(SNR));
    VI_P_CA_B=zeros(1,length(SNR));
    VI_P_SO_AB=zeros(1,length(SNR));
    detect_rate_LSTMOS=zeros(1,length(SNR));
    LSTMOSVI_P_CA_AB=zeros(1,length(SNR));
    LSTMOSVI_P_GO_AB=zeros(1,length(SNR));
    LSTMOSVI_P_CA_A=zeros(1,length(SNR));
    LSTMOSVI_P_CA_B=zeros(1,length(SNR));
    LSTMOSVI_P_LSTMOS_AB=zeros(1,length(SNR));
    detect_rate_LSTMOSVI=zeros(1,length(SNR));
    %%%%%%%%%%%%%%%%%%%%%%Calculate OS-CFAR threshold factor%%%%%%%%%%%%%%%%%
    syms T
    R = 2*N_reference;
    k = round(R*3/4);
    Pfa_OS = k*nchoosek(R,k)*gamma(T+R-k+1)*gamma(k)/gamma(R+T+1);
    T1_OS = vpasolve(Pfa_OS == P_FA, T);
    T2_OS = double(T1_OS);
    T_OS = T2_OS(T2_OS == abs(T2_OS));
    %%%%%%%%%%%%%%%%%%%%%%ideal%%%%%%%%%%%%%%%%%
    tic
    for jj=1:length(SNR)   
        num_detect=0;
        for mont=1:M_num   
            data_pingfang = getexpdata(N_reference,SNR(jj),noise);
            threshold_CFAR=-log(P_FA);  
            if data_pingfang(ii)>=threshold_CFAR
                num_detect=num_detect+1;
            end
        end
        ideal(jj)=num_detect/M_num;   
        save(['noise' num2str(noise) '/ideal.mat'],'ideal');
    end
    disp('Completed ideal');
    elapsedTime=elapsedTime+toc;
    disp(['time:', num2str(elapsedTime/60), ' min']);
    %%%%%%%%%%%%%%%%%%%%%%CA-CFAR%%%%%%%%%%%%%%%%%
    tic
    for jj=1:length(SNR)   
        num_detect=0;
        for mont=1:M_num   
            data_pingfang = getexpdata(N_reference,SNR(jj),noise);
            data_select=sum(data_pingfang(ii-N_reference:ii-1))+sum(data_pingfang(ii+1:ii+N_reference));
            threshold_CFAR=data_select*T_CA;  
            if data_pingfang(ii)>=threshold_CFAR
                num_detect=num_detect+1;
            end
        end
        detect_rate_CA(jj)=num_detect/M_num;   
        save(['noise' num2str(noise) '/detect_rate_CA.mat'],'detect_rate_CA');
    end
    disp('Completed CA-CFAR');
    elapsedTime=elapsedTime+toc;
    disp(['time:', num2str(elapsedTime/60), ' min']);
    %%%%%%%%%%%%%%%%%%%%%%OS-CFAR%%%%%%%%%%%%%%%%%
    tic
    for jj=1:length(SNR)   
        num_detect=0;
        for mont=1:M_num   
            data_pingfang = getexpdata(N_reference,SNR(jj),noise);
            data_select_temp=sort([data_pingfang(ii-N_reference:ii-1) data_pingfang(ii+1:ii+N_reference)]);
            data_select=data_select_temp(k);
            threshold_CFAR=data_select*T_OS;  
            if data_pingfang(ii)>=threshold_CFAR
                num_detect=num_detect+1;
            end
        end
        detect_rate_OS(jj)=num_detect/M_num;   
        save(['noise' num2str(noise) '/detect_rate_OS.mat'],'detect_rate_OS');
    end
    disp('Completed OS-CFAR');
    elapsedTime=elapsedTime+toc;
    disp(['time:', num2str(elapsedTime/60), ' min']);
    %%%%%%%%%%%%%%%%%%%%%%VI-CFAR%%%%%%%%%%%%%%%%%
    tic
    for jj=1:length(SNR)   
        C_CA_AB=0;
        C_GO_AB=0;
        C_CA_A=0;
        C_CA_B=0;
        C_SO_AB=0;
        num_detect=0;
        for mont=1:M_num   
            data_pingfang = getexpdata(N_reference,SNR(jj),noise);
            XA=data_pingfang(ii-N_reference:ii-1);
            XB=data_pingfang(ii+1:ii+N_reference);
            VIA=1+(1/(N_reference-1))*(sum((XA-mean(XA)).^2)/(mean(XA).^2));
            VIB=1+(1/(N_reference-1))*(sum((XB-mean(XB)).^2)/(mean(XB).^2));
            MR=mean(XA)/mean(XB);
            if VIA<=KVI && VIB<=KVI
                if MR>=KMR^-1 && MR<=KMR
                    data_select=sum(data_pingfang(ii-N_reference:ii-1))+sum(data_pingfang(ii+1:ii+N_reference));
                    threshold_CFAR=data_select*T_CA;
                    C_CA_AB=C_CA_AB+1;
                else
                    data_select=max(sum(data_pingfang(ii-N_reference:ii-1)),sum(data_pingfang(ii+1:ii+N_reference)));
                    threshold_CFAR=data_select*T_GOCA;
                    C_GO_AB=C_GO_AB+1;
                end
            elseif VIA>KVI && VIB<=KVI
                data_select=sum(data_pingfang(ii+1:ii+N_reference));
                threshold_CFAR=data_select*T_SOCA;
                C_CA_B=C_CA_B+1;
            elseif VIA<=KVI && VIB>KVI
                data_select=sum(data_pingfang(ii-N_reference:ii-1));
                threshold_CFAR=data_select*T_SOCA;
                C_CA_A=C_CA_A+1;
            else
                data_select=min(sum(data_pingfang(ii-N_reference:ii-1)),sum(data_pingfang(ii+1:ii+N_reference)));
                threshold_CFAR=data_select*T_SOCA;
                C_SO_AB=C_SO_AB+1;
            end
            if data_pingfang(ii)>=threshold_CFAR
                num_detect=num_detect+1;
            end
        end
        VI_P_CA_AB(jj)=C_CA_AB/M_num;
        VI_P_GO_AB(jj)=C_GO_AB/M_num;
        VI_P_CA_A(jj)=C_CA_A/M_num;
        VI_P_CA_B(jj)=C_CA_B/M_num;
        VI_P_SO_AB(jj)=C_SO_AB/M_num;
        detect_rate_VI(jj)=num_detect/M_num;   
        save(['noise' num2str(noise) '/P_VI.mat'],'VI_P_CA_AB','VI_P_GO_AB','VI_P_CA_A','VI_P_CA_B','VI_P_SO_AB');
        save(['noise' num2str(noise) '/detect_rate_VI.mat'],'detect_rate_VI');
    end
    disp('Completed VI-CFAR');
    elapsedTime=elapsedTime+toc;
    disp(['time:', num2str(elapsedTime/60), ' min']);
    %%%%%%%%%%%%%%%%%%%%%%COS-LSTM-CFAR%%%%%%%%%%%%%%%%%
    for jj=1:length(SNR)   
        tic
        num_detect=0;
        for mont=1:M_num   
            data_pingfang = getexpdata(N_reference,SNR(jj),noise);

            X1=[data_pingfang(ii-N_reference:ii-1)];
            minX1=min(X1);
            maxX1=max(X1);
            X1n=(X1-minX1)/(maxX1-minX1)*1000;
            Y1=classify(net1.net, X1n');
            X2=[data_pingfang(ii+1:ii+N_reference)];
            minX2=min(X2);
            maxX2=max(X2);
            X2n=(X2-minX2)/(maxX2-minX2)*1000;
            Y2=classify(net1.net, X2n');
            
            K1=18;
            K2=18;
            while Y1==categorical(2)
                [SX1,I1]=sort(X1);
                X1(I1(18))=min(X1);
                minX1=min(X1);
                maxX1=max(X1);
                if minX1==maxX1
                    break
                end
                X1n=(X1-minX1)/(maxX1-minX1)*1000;
                Y1=classify(net1.net, X1n');
                K1=K1-1;
            end
            while Y2==categorical(2)
                [SX2,I2]=sort(X2);
                X2(I2(18))=min(X2);
                minX2=min(X2);
                maxX2=max(X2);
                if minX2==maxX2
                    break
                end
                X2n=(X2-minX2)/(maxX2-minX2)*1000;
                Y2=classify(net1.net, X2n');
                K2=K2-1;
            end

            syms T
            R = K1+K2;                  
            k = floor(R*3/4);                  
            Pfa_OS = k*nchoosek(R,k)*gamma(T+R-k+1)*gamma(k)/gamma(R+T+1);          
            T1_OS = vpasolve(Pfa_OS == P_FA, T);       
            T2_OS = double(T1_OS);
            T_OS = T2_OS(T2_OS == abs(T2_OS));     

            data_select_temp=sort([data_pingfang(ii-N_reference:ii-1) data_pingfang(ii+1:ii+N_reference)]);  
            data_select=data_select_temp(k);
            threshold_CFAR=data_select*T_OS;  
            if data_pingfang(ii)>=threshold_CFAR
                num_detect=num_detect+1;
            end
            if mod(mont,1000)==0
                disp(['Completed LSTMOS-CFAR：' num2str(jj) '(' num2str(mont) ')' '/' num2str(length(SNR))]);
            end
        end
        detect_rate_LSTMOS(jj)=num_detect/M_num;   
        save(['noise' num2str(noise) '/detect_rate_LSTMOS.mat'],'detect_rate_LSTMOS');
        disp(['Completed LSTMOS-CFAR：' num2str(jj) '/' num2str(length(SNR))]);
        elapsedTime=elapsedTime+toc;
        disp(['time:', num2str(elapsedTime/60), ' min']);
    end
    %%%%%%%%%%%%%%%%%%%%%%ACA-LSTM-CFAR%%%%%%%%%%%%%%%%%
    for jj=1:length(SNR)   
        tic
        C_CA_AB=0;
        C_GO_AB=0;
        C_CA_A=0;
        C_CA_B=0;
        C_LSTMOS_AB=0;
        num_detect=0;
        for mont=1:M_num   
            data_pingfang = getexpdata(N_reference,SNR(jj),noise);
            X1=[data_pingfang(ii-N_reference:ii-1)];
            minX1=min(X1);
            maxX1=max(X1);
            X1n=(X1-minX1)/(maxX1-minX1)*1000;
            Y1=classify(net1.net, X1n');
            X2=[data_pingfang(ii+1:ii+N_reference)];
            minX2=min(X2);
            maxX2=max(X2);
            X2n=(X2-minX2)/(maxX2-minX2)*1000;
            Y2=classify(net1.net, X2n');
            X=[X1 X2];
            minX=min(X);
            maxX=max(X);
            Xn=(X-minX)/(maxX-minX)*1000;
            XA=data_pingfang(ii-N_reference:ii-1);
            XB=data_pingfang(ii+1:ii+N_reference);
            VIA=1+(1/(N_reference-1))*(sum((XA-mean(XA)).^2)/(mean(XA).^2));
            VIB=1+(1/(N_reference-1))*(sum((XB-mean(XB)).^2)/(mean(XB).^2));
            MR=mean(X1)/mean(X2);
            if Y1==categorical(1)&&Y2==categorical(1)
                MRnet=classify(net2.net, Xn');
                if MRnet==categorical(1)
                    data_select=sum(data_pingfang(ii-N_reference:ii-1))+sum(data_pingfang(ii+1:ii+N_reference));
                    threshold_CFAR=data_select*T_CA;
                    C_CA_AB=C_CA_AB+1;
                else
                    data_select=max(sum(data_pingfang(ii-N_reference:ii-1)),sum(data_pingfang(ii+1:ii+N_reference)));
                    threshold_CFAR=data_select*T_GOCA;
                    C_GO_AB=C_GO_AB+1;
                end
            elseif Y1==categorical(2)&&Y2==categorical(1)
                data_select=sum(data_pingfang(ii+1:ii+N_reference));
                threshold_CFAR=data_select*T_SOCA;
                C_CA_B=C_CA_B+1;
            elseif Y1==categorical(1)&&Y2==categorical(2)
                data_select=sum(data_pingfang(ii-N_reference:ii-1));
                threshold_CFAR=data_select*T_SOCA;
                C_CA_A=C_CA_A+1;
            else
                K1=18;
                K2=18;
                while Y1==categorical(2)
                    [SX1,I1]=sort(X1);
                    X1(I1(18))=min(X1);
                    minX1=min(X1);
                    maxX1=max(X1);
                    if minX1==maxX1
                        break
                    end
                    X1n=(X1-minX1)/(maxX1-minX1)*1000;
                    Y1=classify(net1.net, X1n');
                    K1=K1-1;
                end
                while Y2==categorical(2)
                    [SX2,I2]=sort(X2);
                    X2(I2(18))=min(X2);
                    minX2=min(X2);
                    maxX2=max(X2);
                    if minX2==maxX2
                        break
                    end
                    X2n=(X2-minX2)/(maxX2-minX2)*1000;
                    Y2=classify(net1.net, X2n');
                    K2=K2-1;
                end
                syms LSTM_T
                LSTM_R = K1+K2;                 
                K = floor(LSTM_R*3/4);                 
                Pfa_LSTMOS = K*nchoosek(LSTM_R,K)*gamma(LSTM_T+LSTM_R-K+1)*gamma(K)/gamma(LSTM_R+LSTM_T+1);           
                T1_LSTMOS = vpasolve(Pfa_LSTMOS == P_FA, LSTM_T);       
                T2_LSTMOS = double(T1_LSTMOS);
                T_LSTMOS = T2_LSTMOS(T2_LSTMOS == abs(T2_LSTMOS));      

                data_select_temp=sort([data_pingfang(ii-N_reference:ii-1) data_pingfang(ii+1:ii+N_reference)]); 
                data_select=data_select_temp(K);
                threshold_CFAR=data_select*T_LSTMOS;  
                C_LSTMOS_AB=C_LSTMOS_AB+1;
            end
            if data_pingfang(ii)>=threshold_CFAR
                num_detect=num_detect+1;
            end
            if mod(mont,1000)==0
                disp(['Completed LSTMOSVI：' num2str(jj) '(' num2str(mont) ')' '/' num2str(length(SNR))]);
            end
        end
        LSTMOSVI_P_CA_AB(jj)=C_CA_AB/M_num;
        LSTMOSVI_P_GO_AB(jj)=C_GO_AB/M_num;
        LSTMOSVI_P_CA_A(jj)=C_CA_A/M_num;
        LSTMOSVI_P_CA_B(jj)=C_CA_B/M_num;
        LSTMOSVI_P_LSTMOS_AB(jj)=C_LSTMOS_AB/M_num;
        detect_rate_LSTMOSVI(jj)=num_detect/M_num;   
        save(['noise' num2str(noise) '/P_LSTMOSVI.mat'],'LSTMOSVI_P_CA_AB','LSTMOSVI_P_GO_AB','LSTMOSVI_P_CA_A','LSTMOSVI_P_CA_B','LSTMOSVI_P_LSTMOS_AB');
        save(['noise' num2str(noise) '/detect_rate_LSTMOSVI.mat'],'detect_rate_LSTMOSVI');
        disp(['Completed ACA:' num2str(jj) '/' num2str(length(SNR))]);
        elapsedTime=elapsedTime+toc;
        disp(['time:', num2str(elapsedTime/60), ' min']);
    end
end
%% plot
SNR_dB=0:40;
figure(1)
plot(SNR_dB,ideal,'-',LineWidth=1,MarkerSize=8);
hold on;
plot(SNR_dB,detect_rate_CA,'-*',LineWidth=1,MarkerSize=8);
hold on;
plot(SNR_dB,detect_rate_OS,'-diamond',LineWidth=1,MarkerSize=8);
hold on;
plot(SNR_dB,detect_rate_VI,'-o',LineWidth=1,MarkerSize=8);
hold on;
plot(SNR_dB,detect_rate_LSTMOS,'-<',LineWidth=1,MarkerSize=8);
hold on;
plot(SNR_dB,detect_rate_LSTMOSVI,'-pentagram',LineWidth=1,MarkerSize=8);
hold off;
xlabel('\fontname{Times new roman}SNR(dB)',FontSize=15);
ylabel('\fontname{Times new roman}Probability of detection',FontSize=15);
grid on;
legend("ideal","CA","OS","VI","COS-LSTM","ACA-LSTM",fontsize=15,FontName='Times New Roman')
ax = gca;
set(ax, FontName='Times New Roman',fontsize=15);
xlim([0,30])
ylim([0,1])
xticks(0:5:30)
yticks(0:0.1:1)
