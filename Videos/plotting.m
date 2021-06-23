figure(5)
aviobj = avifile('fft.avi', 'fps', 6);
for n= 1: plotPoints + 1
    h = plot(time, InTime(:, n));
    ylim([0 1]);
    set(h,'EraseMode','xor');
    frame = getframe(gca);
    aviobj = addframe(aviobj,frame);
end
aviobj = close(aviobj);