%{ 
File: svd_images.m
Author: Owen Morehead
Purpose: Reads and plots all compressed images from .dat files outputted
form Fortran program as well as the error between compressed and original
image.
%}

%-------- Using data from Fortran to plot SVD compressed images -----

img = importdata('dog_bw_data.dat'); %data from online image website

ks = ["00020" "00040" "00080" "00160" "00320" "00640" "01280" "02560" "03355"];
sigmas = [20,40,80,160,320,640,1280,2560,3355];
%ks = ["00040" "00080" "00320" "01280"];
errors = [];


figure(1); hold on;
sgtitle("Image Compression Using k Singular Values",'fontsize',20,'interpreter','latex')

for i = 1:length(ks)
    
file_name = strcat('Image_appn_1',ks(i),'.dat');
disp(file_name)
img_i = importdata(file_name);

img_i = im2double(img_i); %convert to double for error calculations
errors(i) = norm(img - img_i,'fro')/(size(img,1)*size(img,2));

img_i = uint8(img_i); %convert double to unsigned 8-bit data type for plotting

subplot(3,3,i)
imshow(img_i)
%imagesc(img_i)
%daspect([2 2 2])
%colormap(gca,gray)
set(gca,'XTick',[]); set(gca,'YTick',[])
title(sprintf("k = %.0f",sigmas(i)),'interpreter','latex')
 
    
end


%---can use error data from fortran, svd_errors, or error values calcualted in
%matlab, errors ---

svd_errors = importdata("svd_errors.dat"); 

figure(2); hold on; grid on;
%plot error as a function of number of singular values
plot(sigmas,svd_errors)
xlabel('Number of Singular Values Used', 'Fontsize', 18, 'Interpreter', 'latex')
ylabel('Error Between Compressed and Original Image', 'Fontsize', 18, 'Interpreter', 'latex')
title('Compressed Image Errors','Fontsize', 20, 'Interpreter', 'latex')






