function demo_fig1
    
    cmap =brewermapNL(500, '-RdBu',1.2);
    mygray = [1,1,1] .* 0.8;
    myblue = cmap(50,:);
    myred = cmap(430,:);
    mygreen = [0, 160, 62] ./ 255 .* 1.5;



    dx = 0.01;
    x = (0:dx:10);
    
    s1 = 5 .* exp(-x ./ 4);
    
    s2 = 3 .* exp(-x ./ 4);
    
    
    n1 = x .* 0 + 2;
    n2 = x .* 0 + 1;
    
    [~, x1] = min(abs(s1 - n1));
    x1 = x1 .* dx;
    
    [~, x2] = min(abs(s1 - n2));
    x2 = x2 .* dx;
    
   [~, x3] = min(abs(s2 - n1));
    x3 = x3 .* dx;
 
      [~, x4] = min(abs(s2 - n2));
    x4 = x4 .* dx;
 
    
    
     fig = figure;
    
    fill([10,0,x],[0,0,min([s1;n1])],mygreen)
    hold on
  %  fill([10,0,x],[0,0,min([s1;n2])],'y')

  plot(x1,2, 'ko','MarkerFaceColor','black','MarkerSize',3)
    
    plot(x,n1,'Color',myblue,'LineWidth',1)
    plot(x,n2,'--','Color',myblue,'LineWidth',1)
    plot(x,s1,'Color',myred,'LineWidth',1)
    plot(x,s2,'--','Color',myred,'LineWidth',1)
    
    plot([x1 x1],[0 3],'k:')
  %  plot([x2 x2],[0 5],'k:')
    xlabel("coefficient number n")
    ylabel("contribution \chi_n^2 to MSE")
    
    text(7.1,1.8,"\sigma^2, high noise")
    text(7.1,1.2,"\sigma^2, low noise")
    
    text(1.2,4,"|c_n|^2, bright molecule")
    text(0.3,3,"dimmer")
    
    %text(x1,0.2,'m_{opt}')
    %text(x2,0.2,'later')
    
    xticks([0, x3, x1,x4, x2])
    xticklabels({'0','A','m_{opt}','B','C'})
    yticks([0])
    
        scifig( ...
            "width", 8.2, ...
            "height", 5, ...
            "font", "Arial", ...
            "fontsize", 7, ...
            "tickdir", "in", ...
            "linewidth", 0.5);
    
    
        %    Export the figure
       exportgraphics(fig, "fig_demo.pdf", ...
           "ContentType", "vector"); 
       
     exportgraphics(fig, "fig_demo.png"); 
end