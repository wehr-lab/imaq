vidin = VideoReader('M-0.mj2');
vid = read(vidin);
vid2 = squeeze(vid);

for i=1:size(vid2,3)
    vid2(1:512,1:512,i) = vid2(1:512,1:2:1024,i);
    %NewM(:,:,:,i) = im2uint8(NewM(:,:,:,i));
end
vid2 = vid2(1:512,1:512,:);



%[frameT out] = readTifStandalone;
m = mean(vid2,3);
v = VideoWriter('M-0_dfof','Uncompressed AVI');
v.FrameRate=10;
v.Height = 512;
v.Width=512;
open(v);
figure
colormap(bone(512))
for i = 1:size(vid2,3)
    dfof(1:512,1:512) = 512*(double(vid2(:,:,i))-m)./m+256;
    %dfof=4096*(frame-200)/1.3e+03;
    image(dfof)
    mov = getframe;
    writeVideo(v,mov);
    if mod(i,20)==0
        i
    end
end
close(v);
clear m

vidin = VideoReader('M-0_dfof.avi');
vout2 = read(vidin);
vout3 = (vout2-100)*(256/30);
for j=1:18
    for i=1:20
        avgV(:,:,1,((j-1)*20)+i)=mean(vout3(:,:,1,((j-1)*20)+i:360:end), 4);
        if i<=5
            avgV(1:10,1:10,1,((j-1)*20)+i)=256;
        end
    end
    j
end
v = VideoWriter('M-0_dfof_avg2','Uncompressed AVI');
v.FrameRate=10;
open(v);
% writeVideo(v,uint8(avgV));
% close(v);
figure
    colormap(bone(256));
for i = 1:size(avgV,4);
    image(avgV(:,:,1,i))
    if mod(i,60)<20
        text(12,8, sprintf('%d''', (ceil(i/60))), 'color', 'w', 'fontsi', 18)
    elseif mod(i,60)<40
        text(12,8, sprintf('%d''''', (ceil(i/60))), 'color', 'w', 'fontsi', 18)
    elseif mod(i,60)<60
        text(12,8, sprintf('%d''''''', (ceil(i/60))), 'color', 'w', 'fontsi', 18)
    end
    axis equal;
    axis off
    mov = getframe(gca);
    writeVideo(v,mov);
end
close(v);

vid = VideoWriter('0831_binned');
vid.FrameRate=50;
open(vid);
writeVideo(vid,mov);
close(vid)