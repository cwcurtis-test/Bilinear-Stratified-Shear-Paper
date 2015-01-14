function Movie_Maker_1d(movie_plot,X,no_of_plots,file_name)

title = strcat(file_name,'.avi');

writerObj = VideoWriter(title);
open(writerObj);

plot(X,squeeze(movie_plot(1,:,1))+1,'black',X,squeeze(movie_plot(1,:,2)),'black')
axis([min(X) max(X) -1 2])

set(gca,'nextplot','replacechildren');

for mm = 1:no_of_plots
    
    plot(X,squeeze(movie_plot(mm,:,1))+1,'black',X,squeeze(movie_plot(mm,:,2)),'black')
    axis([min(X) max(X) -1 2])

    frame = getframe;
    writeVideo(writerObj,frame);

end

close(writerObj);
