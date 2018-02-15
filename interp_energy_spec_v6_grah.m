function a=interp_energy_spec_v6
%%indepedent program  
%%  average with time 18 to 23
f_name=['/Volumes/Untitled/1_km/wrfout_d01_0001-01-01_00:00:00';'/Volumes/Untitled/500m/wrfout_d01_0001-01-01_00:00:00';'/Volumes/Untitled/250m/wrfout_d01_0001-01-01_00:00:00'];

dim=[[576,144,101,1];[1152,288,101,1];[2304,576,101,1]]
time=[18 19 20 21 22 23];

sizedim=size(dim(:,1));
%sizedim=1
for ds=1:sizedim(1)
y=[0:dim(ds,2)]
y=0;
for dt=1:6

w_in=ncread(f_name(ds,:),'W',[1,1,1,time(dt)],dim(ds,:)); %% store as string matrix
ph=ncread(f_name(ds,:),'PH',[1,1,1,time(dt)],dim(ds,:));
phb=ncread(f_name(ds,:),'PHB',[1,1,1,time(dt)],dim(ds,:));
z=(ph+phb)/9.8;

h=5000;  %%%height
[nx,ny,nz,nt]=size(w_in);
w=w_in;
sik=size(w)
%w=w_in(:,:,:,time);
%z=z(:,:,:,time);
p=0;
W_XY=zeros(size(z));
W_XY=W_XY(:,:,1);
%%%%%%%%%%%%%%%%%%%interpolation W_XY at height h
for i=1:nx
   for j=1:ny
       p=0;
       for k=1:nz-3
     if (z(i,j,k)<=h && z(i,j,k+1)>=h)
         p=(h-z(i,j,k))/(z(i,j,k+1)-z(i,j,k));
         W_XY(i,j)=w(i,j,k)+(w(i,j,k+1)-w(i,j,k))*p;
         k;
      end  
       end
   end
end





%[nx,ny,nz]=size(wq);
avgy_wq=zeros(nx);
flu_wq=zeros(nx,ny);
sumy_wq=zeros(nx);

%%%%%%%%%%%%%%points at 5km
%pt = find(zq>5000 & zq<5025);
%k_index=floor(pt/nx/ny)+1;
%j_index=floor((pt-(k_index-1)*nx*ny)/ny)+1;
%i_index=pt-(k_index-1)*nx*ny-(j_index-1)*ny;



for i=1:nx
    avgy_wq(i)=sum(W_XY(i,:))/ny;
end 

for j=1:ny
    for i=1:nx
   %flu_w(:,j,:)=w(:,j,:)-avgy_w
   flu_wq(i,j)=W_XY(i,j)-avgy_wq(i);
    end
end

%%%%%%%%%%%%%finding the max flu point.
flu_w_line=flu_wq.^2;

for i=1:nx
sumy_wq(i)=sum(flu_w_line(i,:));
end

max_sumy_wq=0;
for i=1:nx    
      if sumy_wq(i)>max_sumy_wq
          maxi=i
          max_sumy_wq=sumy_wq(i);
      end
end
%maxi=find(sumy_wq==max(sumy_wq))


%%%%%%%%%%%%%%%%%%%%%%%FFT
uy1=W_XY(maxi,:);
N=size(uy1);
N=N(2);

uhat1 = fft(uy1)/N;
x=sum(uy1.^2)/N-sum(abs(uhat1).^2);
%x=sum(uy1.^2)-sum(abs(uhat1).^2)/N;
Ek1 = 0.5*abs(uhat1(1:N/2+1)).^2;
Ek1(2:N/2) = Ek1(2:N/2) + 0.5*abs(uhat1(N:-1:N/2+2)).^2;
y=y+Ek1/6;



%k=[0:N/2];
%loglog(k,a)
%logkhalf,b,'b',logkhalf,c,'g');
%title('energy spectrum');
%xlabel('k')
%ylabel('Ek_w')
end %%%end of time loop
k=[0:dim(ds,2)/2];

loglog(k,y)

%plot(log(k),log(y))

title('energy spectrum');
xlabel('k')
ylabel('Ek_w')
hold on
%clear y
end %%end of dataset loop
hold off


end
