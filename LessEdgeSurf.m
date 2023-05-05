function [] = LessEdgeSurf(surf_plot,xnumlines,ynumlines)
% reduces the amount of mesh lines in a surf plot.
% Input: surf_plot - surface handle.

% turn sample-based edges off:
surf_plot.EdgeColor = 'none';
%%Extract X,Y and Z data from surface plot
X=surf_plot.XData;
Y=surf_plot.YData;
Z=surf_plot.ZData;

%%Divide the lengths by the number of lines needed
xlength = size(Z,2);
ylength = size(Z,1);
xspacing = round(xlength/xnumlines);
yspacing = round(ylength/ynumlines);

%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:ylength
    if i == 1
%         mesh([X(i,:);X(i,:)], [Y(i,:);Y(i,:)], [Z(i,:);Z(i,:)],'EdgeColor',[0.4940 0.1840 0.5560],'Linewidth',1.5);
        mesh([X(i,:);X(i,:)], [Y(i,:);Y(i,:)], [Z(i,:);Z(i,:)],'EdgeColor','green','Linewidth',2.5,...
            'FaceColor','none');
    else
        mesh([X(i,:);X(i,:)], [Y(i,:);Y(i,:)], [Z(i,:);Z(i,:)],'EdgeColor',[0 0 0],...
            'FaceColor','none');
    end
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:xlength
  if i == 1
      mesh([X(:,i),X(:,i)], [Y(:,i),Y(:,i)], [Z(:,i),Z(:,i)],'EdgeColor','red','Linewidth',2.5);
  else
      mesh([X(:,i),X(:,i)], [Y(:,i),Y(:,i)], [Z(:,i),Z(:,i)],'EdgeColor',[0 0 0]);
  end
end
hold off
end

