points = [[0.0,-0.5]; [0.5,-0.5]; [0.5,0]; [0.46194,-0.191342]; ...
  [0.353553,-0.353553]; [0.191342,-0.46194]];
mesh = [[2 0 4]; [5 4 0]; [5 0 1]; [5 1 4]; [3 2 4]; [3 4 1]; [3 1 2]];
mesh = mesh+1;
figure, hold on
plot(points(:,1),points(:,2),'or')
for k=1:length(mesh)
  tmpx = points(mesh(k,:),1); tmpx(end+1) = tmpx(1);
  tmpy = points(mesh(k,:),2); tmpy(end+1) = tmpy(1);
  plot(tmpx, tmpy)
end
hold off