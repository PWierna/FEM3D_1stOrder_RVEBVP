ncon_21 = MODEL.Conectivity(21,3:end);
plot3(MODEL.Coordinates(ncon_21,2),MODEL.Coordinates(ncon_21,3),MODEL.Coordinates(ncon_21,4),'ok','LineWidth',1.5);hold on;grid on
xlabel('X');ylabel('Y');zlabel('Z')

for n=1:27
    text(MODEL.Coordinates(ncon_21(n),2),MODEL.Coordinates(ncon_21(n),3),MODEL.Coordinates(ncon_21(n),4),['\leftarrow ',char(string(n))],'FontSize',15)
end

%plot face elements:
ef_2Dnodes  = MODEL.Conectivity( e , 2 + faces_connect(ef,3:end) );
ef_2Dnodes  = ef_2Dnodes( faceidxs_resort2D(lface_idx,:) );
plot3(MODEL.Coordinates(ef_2Dnodes,2),MODEL.Coordinates(ef_2Dnodes,3),MODEL.Coordinates(ef_2Dnodes,4),'or','LineWidth',1.5);hold on;grid on

for nf=1:length(ef_2Dnodes)
    text(MODEL.Coordinates(ef_2Dnodes(nf),2),MODEL.Coordinates(ef_2Dnodes(nf),3),MODEL.Coordinates(ef_2Dnodes(nf),4),['\uparrow ',char(string(nf))],'FontSize',15,'Color','r')
end

figure(2);
for nf=1:length(ef_2Dnodes)
    plot(ef_2Dcoords(nf,1),ef_2Dcoords(nf,2),'or','LineWidth',1.5);hold on
    text(MODEL.Coordinates(ef_2Dnodes(nf),2),MODEL.Coordinates(ef_2Dnodes(nf),3),MODEL.Coordinates(ef_2Dnodes(nf),4),['\uparrow ',char(string(nf))],'FontSize',15,'Color','r')
end
xlabel('X');ylabel('Y');zlabel('Z')