% question.1.4.3
% detect people on a crosswalk video

clc; clear;
v=VideoReader('crossWalkSelectedFrames.avi') % just 30 frames of the orginal file
t=0;

while hasFrame(v)
    filename = 'saved_gif.gif';
    I = readFrame(v);
    del = 0.001;
    % Your code for person detection
    imshow(detect(I));
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if t == 0
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);         
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
    end
    t=t+1;
end


function I = detect(frame)

    thresh = 1.7; % example
    peopleDetector = vision.PeopleDetector('UprightPeople_96x48','ClassificationThreshold',thresh);
    I = frame;
    [bboxes,scores] = peopleDetector(I);
    I = insertObjectAnnotation(I,'rectangle',bboxes,scores);
end