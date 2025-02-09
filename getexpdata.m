function data_pingfang = getexpdata(N_reference,SNR,noise)
%Obtain exponential distribution echoes
noise_P=1;
data_pingfang=[exprnd(noise_P,1,N_reference) exprnd(noise_P).*SNR exprnd(noise_P,1,N_reference)];
switch noise %Add interference
    case 0

    case 1
        data_pingfang((N_reference+1)-9)=exprnd(noise_P).*SNR;
    case 2
        data_pingfang((N_reference+1)-9)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-4)=exprnd(noise_P).*SNR;
    case 3
        data_pingfang((N_reference+1)+9)=exprnd(noise_P).*SNR;
    case 4
        data_pingfang((N_reference+1)+4)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+9)=exprnd(noise_P).*SNR;
    case 5
        data_pingfang((N_reference+1)-9)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+9)=exprnd(noise_P).*SNR;
    case 6
        data_pingfang((N_reference+1)-9)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-4)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+4)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+9)=exprnd(noise_P).*SNR;
    case 7
        data_pingfang((N_reference+1)-16)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-12)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-8)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-4)=exprnd(noise_P).*SNR;
    case 8
        data_pingfang((N_reference+1)+4)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+8)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+12)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+16)=exprnd(noise_P).*SNR;
    case 9
        data_pingfang((N_reference+1)-16)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-12)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-8)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)-4)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+4)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+8)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+12)=exprnd(noise_P).*SNR;
        data_pingfang((N_reference+1)+16)=exprnd(noise_P).*SNR;
end
end

