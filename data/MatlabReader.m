figure('Renderer', 'painters', 'Position', [10 10 1400 600])

% Take the input data from SimulationParameters.txt
Lx = 40; % d_i
Ly = 20; % d_i

nxc = 128
nyc = 64

x = linspace(0,Lx,nxc)
y = linspace(0,Ly,nyc)


for i=10:10:3000
    index_string = sprintf('%d',i)
    
    %colormap(redblue)
    colormap(turbo)
    
    subplot(2,3,1)
    Ex_filename = "Ex_" + index_string + ".spic"
    Ex = load(Ex_filename);
    pcolor(x,y,Ex)
    titleEx = "E_x - cycle: " + index_string 
    title(titleEx)
    colorbar
    xlabel('x (d_i)')
    ylabel('y (d_i)')

    subplot(2,3,2)
    Ey_filename = "Ey_" + index_string + ".spic"
    Ey = load(Ey_filename);
    pcolor(x,y,Ey)
    titleEy = "E_y - cycle: " + index_string 
    title(titleEy)
    colorbar
    xlabel('x (d_i)')
    ylabel('y (d_i)')

    
    subplot(2,3,3)
    Ez_filename = "Ez_" + index_string + ".spic"
    Ez = load(Ez_filename);
    pcolor(x,y,Ez)
    titleEz = "E_z - cycle: " + index_string 
    title(titleEz)
    colorbar
    xlabel('x (d_i)')
    ylabel('y (d_i)')



    subplot(2,3,4)
    Bx_filename = "Bx_" + index_string + ".spic"
    Bx = load(Bx_filename);
    pcolor(x,y,Bx)
    titleBx = "B_x - cycle: " + index_string 
    title(titleBx)
    colorbar
    xlabel('x (d_i)')
    ylabel('y (d_i)')

    subplot(2,3,5)
    By_filename = "By_" + index_string + ".spic"
    By = load(By_filename);
    pcolor(x,y,By)
    titleBx = "B_y - cycle: " + index_string 
    title(titleBx)
    colorbar
    xlabel('x (d_i)')
    ylabel('y (d_i)')

    
    subplot(2,3,6)
    Bz_filename = "Bz_" + index_string + ".spic"
    Bz = load(Bz_filename);
    pcolor(x,y,Bz)
    titleBz = "B_z - cycle: " + index_string 
    title(titleBz)
    colorbar
    xlabel('x (d_i)')
    ylabel('y (d_i)')
    
    
    
    
    % pause a little bit :)
    pause(0.1)
end