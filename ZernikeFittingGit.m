%% 基本准备
% 检测波长
lambda = 632.8*1E-9; % 单位：m
% 面形采样间隔
SamplegGap = 0.005;
[x,y] = meshgrid(-1 : SamplegGap : 1);

% 圆域内作图
xcircle = x(x.^2+y.^2<1);
ycircle = y(x.^2+y.^2<1);
% 三角剖分
T = delaunay(xcircle,ycircle);

% 计算Zernike多项式
% 1可以符号计算
% 2也可以手敲^_^
% 存贮Zernike多项式的大元胞数组，项数Num_j
Zernike = cell(Num_j,1);
% Zernike第一项是常值1函数
Zernike{1} = ones(size(x));
% 从第二项开始计算
for i = 2 : Num_j
    % 采用符号计算后，可以将其转换成handle方便带入数值
    Zernike{i} = Zernike_handles{i}(x,y);
    flag = Zernike{i}(x.^2+y.^2==0);
    % 中间那个元值如果是nan，那就赋值为0，因为洛必达法则[符号生成函数的原因，手敲可以忽略]
    if isnan(flag)
        Zernike{i}(x.^2+y.^2==0) = 0;
    end
end

%% 点集正交化
% 面形
V = load('面形.mat');
randZ = V.z;
% 施密特正交化得到点集正交多项式表
PointOthZernikeVec = cell(Num_j,1);

% 向量化有效圆域正交Zernike多项式作为点集正交化初值
% FaceVector[是一个根据自拟规则将有效点向量化的函数]
[X,Y,PointOthZernikeVec{1}] = FaceVectorize(x,y,Zernike{1},1) ;

% 预留点集正交化系数
D = zeros(Num_j);
for j = 2 : Num_j
    % Zernike第j项有效点向量化
    % 1.赋给点集正交多项式做初值
    [~,~,PointOthZernikeVec{j}] = FaceVectorize(x,y,Zernike{j},1);
    % 2.留着求正交化系数D_{jp}
    Temp_Zernike = PointOthZernikeVec{j};
    for p = 1 : j-1
        % 计算点集正交化系数
        D(j,p) = Temp_Zernike' * PointOthZernikeVec{p} / (PointOthZernikeVec{p}' * PointOthZernikeVec{p});
        % 累加得到点集正交多项式
        PointOthZernikeVec{j} = PointOthZernikeVec{j} + D(j,p) * PointOthZernikeVec{p};
    end 
end

% 计算面形数据点集正交系数B_{j}
B = zeros(Num_j,1);
[~,~,MeasuredFace] = FaceVectorize(x,y,randZ,1);
for p = 1 : Num_j
    B(p) = MeasuredFace' * PointOthZernikeVec{p} / (PointOthZernikeVec{p}' * PointOthZernikeVec{p});
end

% 计算点集和圆域转换系数C_{ji}
C = zeros(Num_j); 
for j = 2 : Num_j
    C(j,j) = 1;
    for i = 1 : j - 1 
        C(j,i) = D(j,i);
        for s = 1 : j - i
            C(j,i) =C(j,i) + D(j,j-s)*C(j-s,j);
        end
    end
end

% 计算圆域正交多项式系数
A = zeros(Num_j,1);
A(end) = B(end);
for j = 1 : Num_j - 1
    A(j) = B(j);
    for  i = j +1 : Num_j
        A(j) = A(j) + B(i)*C(i,j);
    end
end


%% 验证拟合是否相似，画个动图
fig = figure;
W = 0;
im = cell(Num_j,1);
for i = 1 : Num_j
    W = A(i)*Zernike{i} + W;
    trisurf(T,xcircle,ycircle,W(x.^2+y.^2<1))
    shading interp
    axis equal
    title(sprintf('加上Zernike第%d项',i))
    xlabel('归一化x方向,单位mm')
    ylabel('归一化y方向,单位mm')
    colorbar
    colormap('jet')
    view(0,90)
    pause(0.1)
    frame = getframe(fig);
    im{i} = frame2im(frame);
end
filename = '.\拟合动画.gif';

for idx = 1: Num_j
    % 制作gif文件，图像必须是index索引图像
    [Frame, map] = rgb2ind(im{idx}, 256);
    if idx == 1
        imwrite(Frame, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.3);
    else
        imwrite(Frame, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.3);
    end
end

%% 计算残差时要回归到0，不然无法计算
W(isnan(W))=0;
W(x.^2+y.^2>1) = 0;
randZ(isnan(randZ))=0;
% 只在圆域内作图
figure('Name',sprintf('Zernike%d项拟合残差图',Num_j))
subplot(1,3,1)
% 测量图
trisurf(T,x(x.^2+y.^2<1),y(x.^2+y.^2<1),randZ(x.^2+y.^2<1))
shading interp
colormap('jet')
title('实际测量图')
view(0,90)
colorbar
axis equal
axis([-1 1 -1 1])


subplot(1,3,2)
% 拟合图
trisurf(T,x(x.^2+y.^2<1),y(x.^2+y.^2<1),W(x.^2+y.^2<1))
shading interp
colormap('jet')
title('Zernike拟合图')
view(0,90)
colorbar
axis equal
axis([-1 1 -1 1])

subplot(1,3,3)
% 拟合残差 = 实测-拟合
Res =  randZ - W;
trisurf(T,x(x.^2+y.^2<1),y(x.^2+y.^2<1),Res(x.^2+y.^2<1))
shading interp
colormap('jet')
title('拟合残差图')
view(0,90)
colorbar
axis equal
axis([-1 1 -1 1])