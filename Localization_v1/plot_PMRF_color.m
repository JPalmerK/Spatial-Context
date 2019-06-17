

x=[];
y=[];
times=[];

 figure(1)
 hold on;

   x=[];
y=[];
times=[];

for(j=1:length(array_struct))

x=localize_struct.hyd(array_struct(j).master).coordinates(2,:);
y=localize_struct.hyd(array_struct(j).master).coordinates(1,:);
times=localize_struct.hyd(array_struct(j).master).rtimes(1,:);

scatter(x, y, 10, times/hyd(1).detection.parm.sample_freq/3600, 'filled');

end


hold on
colorbar
grid on



for(j=1:length(hydrophone_struct))
    hold on;
    plot(hydrophone_struct(j).location(2), hydrophone_struct(j).location(1),'k*');
    
end