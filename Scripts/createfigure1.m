function createfigure1(cdata1, X1, Y1, S1, C1, X2, Y2, C2, X3, Y3, C3, cdata2)
%CREATEFIGURE(cdata1, X1, Y1, S1, C1, X2, Y2, C2, X3, Y3, C3, cdata2)
%  CDATA1:  image cdata
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c
%  X2:  scatter x
%  Y2:  scatter y
%  C2:  scatter c
%  X3:  scatter x
%  Y3:  scatter y
%  C3:  scatter c
%  CDATA2:  image cdata

%  Auto-generated by MATLAB on 21-Aug-2019 15:38:10

% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],...
    'Renderer','painters');

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure1);
hold(subplot1,'on');

% Create image
image([-70.5 -70.499 -70.498 -70.497 -70.496 -70.495 -70.494 -70.493 -70.492 -70.491 -70.49 -70.489 -70.488 -70.487 -70.486 -70.485 -70.484 -70.483 -70.482 -70.481 -70.48 -70.479 -70.478 -70.477 -70.476 -70.475 -70.474 -70.473 -70.472 -70.471 -70.47 -70.469 -70.468 -70.467 -70.466 -70.465 -70.464 -70.463 -70.462 -70.461 -70.46 -70.459 -70.458 -70.457 -70.456 -70.455 -70.454 -70.453 -70.452 -70.451 -70.45 -70.449 -70.448 -70.447 -70.446 -70.445 -70.444 -70.443 -70.442 -70.441 -70.44 -70.439 -70.438 -70.437 -70.436 -70.435 -70.434 -70.433 -70.432 -70.431 -70.43 -70.429 -70.428 -70.427 -70.426 -70.425 -70.424 -70.423 -70.422 -70.421 -70.42 -70.419 -70.418 -70.417 -70.416 -70.415 -70.414 -70.413 -70.412 -70.411 -70.41 -70.409 -70.408 -70.407 -70.406 -70.405 -70.404 -70.403 -70.402 -70.401 -70.4 -70.399 -70.398 -70.397 -70.396 -70.395 -70.394 -70.393 -70.392 -70.391 -70.39 -70.389 -70.388 -70.387 -70.386 -70.385 -70.384 -70.383 -70.382 -70.381 -70.38 -70.379 -70.378 -70.377 -70.376 -70.375 -70.374 -70.373 -70.372 -70.371 -70.37 -70.369 -70.368 -70.367 -70.366 -70.365 -70.364 -70.363 -70.362 -70.361 -70.36 -70.359 -70.358 -70.357 -70.356 -70.355 -70.354 -70.353 -70.352 -70.351 -70.35 -70.349 -70.348 -70.347 -70.346 -70.345 -70.344 -70.343 -70.342 -70.341 -70.34 -70.339 -70.338 -70.337 -70.336 -70.335 -70.334 -70.333 -70.332 -70.331 -70.33 -70.329 -70.328 -70.327 -70.326 -70.325 -70.324 -70.323 -70.322 -70.321 -70.32 -70.319 -70.318 -70.317 -70.316 -70.315 -70.314 -70.313 -70.312 -70.311 -70.31 -70.309 -70.308 -70.307 -70.306 -70.305 -70.304 -70.303 -70.302 -70.301 -70.3 -70.299 -70.298 -70.297 -70.296 -70.295 -70.294 -70.293 -70.292 -70.291 -70.29 -70.289 -70.288 -70.287 -70.286 -70.285 -70.284 -70.283 -70.282 -70.281 -70.28 -70.279 -70.278 -70.277 -70.276 -70.275 -70.274 -70.273 -70.272 -70.271 -70.27 -70.269 -70.268 -70.267 -70.266 -70.265 -70.264 -70.263 -70.262 -70.261 -70.26 -70.259 -70.258 -70.257 -70.256 -70.255 -70.254 -70.253 -70.252 -70.251 -70.25 -70.249 -70.248 -70.247 -70.246 -70.245 -70.244 -70.243 -70.242 -70.241 -70.24 -70.239 -70.238 -70.237 -70.236 -70.235 -70.234 -70.233 -70.232 -70.231 -70.23 -70.229 -70.228 -70.227 -70.226 -70.225 -70.224 -70.223 -70.222 -70.221 -70.22 -70.219 -70.218 -70.217 -70.216 -70.215 -70.214 -70.213 -70.212 -70.211 -70.21 -70.209 -70.208 -70.207 -70.206 -70.205 -70.204 -70.203 -70.202 -70.201 -70.2 -70.199 -70.198 -70.197 -70.196 -70.195 -70.194 -70.193 -70.192 -70.191 -70.19 -70.189 -70.188 -70.187 -70.186 -70.185 -70.184 -70.183 -70.182 -70.181 -70.18 -70.179 -70.178 -70.177 -70.176 -70.175 -70.174 -70.173 -70.172 -70.171 -70.17 -70.169 -70.168 -70.167 -70.166 -70.165 -70.164 -70.163 -70.162 -70.161 -70.16 -70.159 -70.158 -70.157 -70.156 -70.155 -70.154 -70.153 -70.152 -70.151 -70.15 -70.149 -70.148 -70.147 -70.146 -70.145 -70.144 -70.143 -70.142 -70.141 -70.14 -70.139 -70.138 -70.137 -70.136 -70.135 -70.134 -70.133 -70.132 -70.131 -70.13 -70.129 -70.128 -70.127 -70.126 -70.125 -70.124 -70.123 -70.122 -70.121 -70.12 -70.119 -70.118 -70.117 -70.116 -70.115 -70.114 -70.113 -70.112 -70.111 -70.11 -70.109 -70.108 -70.107 -70.106 -70.105 -70.104 -70.103 -70.102 -70.101 -70.1],...
    [42 42.001 42.002 42.003 42.004 42.005 42.006 42.007 42.008 42.009 42.01 42.011 42.012 42.013 42.014 42.015 42.016 42.017 42.018 42.019 42.02 42.021 42.022 42.023 42.024 42.025 42.026 42.027 42.028 42.029 42.03 42.031 42.032 42.033 42.034 42.035 42.036 42.037 42.038 42.039 42.04 42.041 42.042 42.043 42.044 42.045 42.046 42.047 42.048 42.049 42.05 42.051 42.052 42.053 42.054 42.055 42.056 42.057 42.058 42.059 42.06 42.061 42.062 42.063 42.064 42.065 42.066 42.067 42.068 42.069 42.07 42.071 42.072 42.073 42.074 42.075 42.076 42.077 42.078 42.079 42.08 42.081 42.082 42.083 42.084 42.085 42.086 42.087 42.088 42.089 42.09 42.091 42.092 42.093 42.094 42.095 42.096 42.097 42.098 42.099 42.1 42.101 42.102 42.103 42.104 42.105 42.106 42.107 42.108 42.109 42.11 42.111 42.112 42.113 42.114 42.115 42.116 42.117 42.118 42.119 42.12 42.121 42.122 42.123 42.124 42.125 42.126 42.127 42.128 42.129 42.13 42.131 42.132 42.133 42.134 42.135 42.136 42.137 42.138 42.139 42.14 42.141 42.142 42.143 42.144 42.145 42.146 42.147 42.148 42.149 42.15 42.151 42.152 42.153 42.154 42.155 42.156 42.157 42.158 42.159 42.16 42.161 42.162 42.163 42.164 42.165 42.166 42.167 42.168 42.169 42.17 42.171 42.172 42.173 42.174 42.175 42.176 42.177 42.178 42.179 42.18 42.181 42.182 42.183 42.184 42.185 42.186 42.187 42.188 42.189 42.19 42.191 42.192 42.193 42.194 42.195 42.196 42.197 42.198 42.199 42.2 42.201 42.202 42.203 42.204 42.205 42.206 42.207 42.208 42.209 42.21 42.211 42.212 42.213 42.214 42.215 42.216 42.217 42.218 42.219 42.22 42.221 42.222 42.223 42.224 42.225 42.226 42.227 42.228 42.229 42.23 42.231 42.232 42.233 42.234 42.235 42.236 42.237 42.238 42.239 42.24 42.241 42.242 42.243 42.244 42.245 42.246 42.247 42.248 42.249 42.25 42.251 42.252 42.253 42.254 42.255 42.256 42.257 42.258 42.259 42.26 42.261 42.262 42.263 42.264 42.265 42.266 42.267 42.268 42.269 42.27 42.271 42.272 42.273 42.274 42.275 42.276 42.277 42.278 42.279 42.28 42.281 42.282 42.283 42.284 42.285 42.286 42.287 42.288 42.289 42.29 42.291 42.292 42.293 42.294 42.295 42.296 42.297 42.298 42.299 42.3],...
    cdata1,'Parent',subplot1,'CDataMapping','scaled');

% Create scatter
scatter(X1,Y1,S1,C1,'Parent',subplot1,'MarkerFaceColor','flat',...
    'MarkerEdgeColor','none',...
    'Marker','diamond');

% Create scatter
scatter(X2,Y2,S1,C2,'Parent',subplot1,'MarkerFaceColor','flat',...
    'MarkerEdgeColor','none',...
    'Marker','diamond');

% Create scatter
scatter(X3,Y3,S1,C3,'Parent',subplot1,'MarkerFaceColor','flat',...
    'MarkerEdgeColor','none',...
    'Marker','diamond');

% Create ylabel
ylabel('Latitude','FontWeight','bold','FontSize',11,'FontName','sansserif');

% Create xlabel
xlabel('Longitude','FontWeight','bold','FontSize',11,'FontName','sansserif');

% Create title
title('Projected Ambiguity Surface','FontWeight','bold','FontSize',11,...
    'FontName','sansserif');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[-70.5005 -70.0995]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot1,[41.9995 42.3005]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot1,[-1 1]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'FontName','sansserif','FontWeight','bold','Layer','top',...
    'XTick',[-70.5 -70.3 -70.1],'XTickLabel',{'-70.5�','-70.3�','-70.1�'},...
    'YTick',[42 42.1 42.2 42.3],'YTickLabel',{'42.0�','42.1�','42.2�','42.3�'});
% Create colorbar
colorbar('peer',subplot1,'FontWeight','bold','FontName','sansserif');

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure1);
hold(subplot2,'on');

% Create image
image([-70.5 -70.499 -70.498 -70.497 -70.496 -70.495 -70.494 -70.493 -70.492 -70.491 -70.49 -70.489 -70.488 -70.487 -70.486 -70.485 -70.484 -70.483 -70.482 -70.481 -70.48 -70.479 -70.478 -70.477 -70.476 -70.475 -70.474 -70.473 -70.472 -70.471 -70.47 -70.469 -70.468 -70.467 -70.466 -70.465 -70.464 -70.463 -70.462 -70.461 -70.46 -70.459 -70.458 -70.457 -70.456 -70.455 -70.454 -70.453 -70.452 -70.451 -70.45 -70.449 -70.448 -70.447 -70.446 -70.445 -70.444 -70.443 -70.442 -70.441 -70.44 -70.439 -70.438 -70.437 -70.436 -70.435 -70.434 -70.433 -70.432 -70.431 -70.43 -70.429 -70.428 -70.427 -70.426 -70.425 -70.424 -70.423 -70.422 -70.421 -70.42 -70.419 -70.418 -70.417 -70.416 -70.415 -70.414 -70.413 -70.412 -70.411 -70.41 -70.409 -70.408 -70.407 -70.406 -70.405 -70.404 -70.403 -70.402 -70.401 -70.4 -70.399 -70.398 -70.397 -70.396 -70.395 -70.394 -70.393 -70.392 -70.391 -70.39 -70.389 -70.388 -70.387 -70.386 -70.385 -70.384 -70.383 -70.382 -70.381 -70.38 -70.379 -70.378 -70.377 -70.376 -70.375 -70.374 -70.373 -70.372 -70.371 -70.37 -70.369 -70.368 -70.367 -70.366 -70.365 -70.364 -70.363 -70.362 -70.361 -70.36 -70.359 -70.358 -70.357 -70.356 -70.355 -70.354 -70.353 -70.352 -70.351 -70.35 -70.349 -70.348 -70.347 -70.346 -70.345 -70.344 -70.343 -70.342 -70.341 -70.34 -70.339 -70.338 -70.337 -70.336 -70.335 -70.334 -70.333 -70.332 -70.331 -70.33 -70.329 -70.328 -70.327 -70.326 -70.325 -70.324 -70.323 -70.322 -70.321 -70.32 -70.319 -70.318 -70.317 -70.316 -70.315 -70.314 -70.313 -70.312 -70.311 -70.31 -70.309 -70.308 -70.307 -70.306 -70.305 -70.304 -70.303 -70.302 -70.301 -70.3 -70.299 -70.298 -70.297 -70.296 -70.295 -70.294 -70.293 -70.292 -70.291 -70.29 -70.289 -70.288 -70.287 -70.286 -70.285 -70.284 -70.283 -70.282 -70.281 -70.28 -70.279 -70.278 -70.277 -70.276 -70.275 -70.274 -70.273 -70.272 -70.271 -70.27 -70.269 -70.268 -70.267 -70.266 -70.265 -70.264 -70.263 -70.262 -70.261 -70.26 -70.259 -70.258 -70.257 -70.256 -70.255 -70.254 -70.253 -70.252 -70.251 -70.25 -70.249 -70.248 -70.247 -70.246 -70.245 -70.244 -70.243 -70.242 -70.241 -70.24 -70.239 -70.238 -70.237 -70.236 -70.235 -70.234 -70.233 -70.232 -70.231 -70.23 -70.229 -70.228 -70.227 -70.226 -70.225 -70.224 -70.223 -70.222 -70.221 -70.22 -70.219 -70.218 -70.217 -70.216 -70.215 -70.214 -70.213 -70.212 -70.211 -70.21 -70.209 -70.208 -70.207 -70.206 -70.205 -70.204 -70.203 -70.202 -70.201 -70.2 -70.199 -70.198 -70.197 -70.196 -70.195 -70.194 -70.193 -70.192 -70.191 -70.19 -70.189 -70.188 -70.187 -70.186 -70.185 -70.184 -70.183 -70.182 -70.181 -70.18 -70.179 -70.178 -70.177 -70.176 -70.175 -70.174 -70.173 -70.172 -70.171 -70.17 -70.169 -70.168 -70.167 -70.166 -70.165 -70.164 -70.163 -70.162 -70.161 -70.16 -70.159 -70.158 -70.157 -70.156 -70.155 -70.154 -70.153 -70.152 -70.151 -70.15 -70.149 -70.148 -70.147 -70.146 -70.145 -70.144 -70.143 -70.142 -70.141 -70.14 -70.139 -70.138 -70.137 -70.136 -70.135 -70.134 -70.133 -70.132 -70.131 -70.13 -70.129 -70.128 -70.127 -70.126 -70.125 -70.124 -70.123 -70.122 -70.121 -70.12 -70.119 -70.118 -70.117 -70.116 -70.115 -70.114 -70.113 -70.112 -70.111 -70.11 -70.109 -70.108 -70.107 -70.106 -70.105 -70.104 -70.103 -70.102 -70.101 -70.1],...
    [42 42.001 42.002 42.003 42.004 42.005 42.006 42.007 42.008 42.009 42.01 42.011 42.012 42.013 42.014 42.015 42.016 42.017 42.018 42.019 42.02 42.021 42.022 42.023 42.024 42.025 42.026 42.027 42.028 42.029 42.03 42.031 42.032 42.033 42.034 42.035 42.036 42.037 42.038 42.039 42.04 42.041 42.042 42.043 42.044 42.045 42.046 42.047 42.048 42.049 42.05 42.051 42.052 42.053 42.054 42.055 42.056 42.057 42.058 42.059 42.06 42.061 42.062 42.063 42.064 42.065 42.066 42.067 42.068 42.069 42.07 42.071 42.072 42.073 42.074 42.075 42.076 42.077 42.078 42.079 42.08 42.081 42.082 42.083 42.084 42.085 42.086 42.087 42.088 42.089 42.09 42.091 42.092 42.093 42.094 42.095 42.096 42.097 42.098 42.099 42.1 42.101 42.102 42.103 42.104 42.105 42.106 42.107 42.108 42.109 42.11 42.111 42.112 42.113 42.114 42.115 42.116 42.117 42.118 42.119 42.12 42.121 42.122 42.123 42.124 42.125 42.126 42.127 42.128 42.129 42.13 42.131 42.132 42.133 42.134 42.135 42.136 42.137 42.138 42.139 42.14 42.141 42.142 42.143 42.144 42.145 42.146 42.147 42.148 42.149 42.15 42.151 42.152 42.153 42.154 42.155 42.156 42.157 42.158 42.159 42.16 42.161 42.162 42.163 42.164 42.165 42.166 42.167 42.168 42.169 42.17 42.171 42.172 42.173 42.174 42.175 42.176 42.177 42.178 42.179 42.18 42.181 42.182 42.183 42.184 42.185 42.186 42.187 42.188 42.189 42.19 42.191 42.192 42.193 42.194 42.195 42.196 42.197 42.198 42.199 42.2 42.201 42.202 42.203 42.204 42.205 42.206 42.207 42.208 42.209 42.21 42.211 42.212 42.213 42.214 42.215 42.216 42.217 42.218 42.219 42.22 42.221 42.222 42.223 42.224 42.225 42.226 42.227 42.228 42.229 42.23 42.231 42.232 42.233 42.234 42.235 42.236 42.237 42.238 42.239 42.24 42.241 42.242 42.243 42.244 42.245 42.246 42.247 42.248 42.249 42.25 42.251 42.252 42.253 42.254 42.255 42.256 42.257 42.258 42.259 42.26 42.261 42.262 42.263 42.264 42.265 42.266 42.267 42.268 42.269 42.27 42.271 42.272 42.273 42.274 42.275 42.276 42.277 42.278 42.279 42.28 42.281 42.282 42.283 42.284 42.285 42.286 42.287 42.288 42.289 42.29 42.291 42.292 42.293 42.294 42.295 42.296 42.297 42.298 42.299 42.3],...
    cdata2,'Parent',subplot2,'CDataMapping','scaled');

% Create scatter
scatter(X1,Y1,S1,C1,'Parent',subplot2,'MarkerFaceColor','flat',...
    'MarkerEdgeColor','none',...
    'Marker','diamond');

% Create scatter
scatter(X2,Y2,S1,C2,'Parent',subplot2,'MarkerFaceColor','flat',...
    'MarkerEdgeColor','none',...
    'Marker','diamond');

% Create scatter
scatter(X3,Y3,S1,C3,'Parent',subplot2,'MarkerFaceColor','flat',...
    'MarkerEdgeColor','none',...
    'Marker','diamond');

% Create zlabel
zlabel('Likelihood','FontWeight','bold','FontName','sansserif');

% Create ylabel
ylabel('Latitude','FontWeight','bold','FontSize',11,'FontName','sansserif');

% Create xlabel
xlabel('Longitude','FontWeight','bold','FontSize',11,'FontName','sansserif');

% Create title
title('Projected Ambiguity Surface','FontWeight','bold','FontSize',11,...
    'FontName','sansserif');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[-70.5005 -70.0995]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot2,[41.9995 42.3005]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot2,[-1 1]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'Colormap',...
    [0.2422 0.1504 0.6603;0.248424761904762 0.161492380952381 0.696258857142857;0.254228571428571 0.17372380952381 0.730246666666667;0.259719428571429 0.186254285714286 0.763479428571429;0.264965333333333 0.198433904761905 0.796860571428571;0.269463809523809 0.211292380952381 0.82814;0.273148952380952 0.225630857142857 0.855755428571429;0.276133714285714 0.241160761904762 0.879973142857143;0.278462666666667 0.257660380952381 0.900913714285714;0.280008 0.274655428571428 0.918415428571428;0.28093619047619 0.291650476190476 0.933665714285714;0.281221523809524 0.308573523809524 0.947319619047619;0.280828571428571 0.325386857142857 0.959540571428571;0.279652380952381 0.342041714285714 0.970021333333333;0.277229714285714 0.358757523809524 0.978862476190476;0.273548571428571 0.375688571428571 0.985982857142857;0.268366857142857 0.392900571428571 0.991329142857143;0.261016571428571 0.410478285714286 0.994792571428571;0.249220380952381 0.428311809523809 0.997656952380952;0.233741523809524 0.446579809523809 0.998152380952381;0.215780952380952 0.465149523809524 0.995659047619048;0.197305714285714 0.483740571428571 0.989477714285714;0.187024761904762 0.501028761904762 0.982415238095238;0.181119047619048 0.517684571428571 0.974208952380952;0.178642857142857 0.528857142857143 0.968157142857143;0.175254212454212 0.55303663003663 0.949534798534799;0.164206593406593 0.576396703296703 0.931541758241758;0.150320879120879 0.598931868131868 0.915364835164835;0.141102564102564 0.620509523809524 0.901354578754579;0.127862637362637 0.641705494505494 0.890407692307692;0.112295604395604 0.662148351648352 0.877239560439561;0.0931835164835166 0.680978021978022 0.85820989010989;0.059824175824176 0.697862637362637 0.834043956043956;0.0196300366300368 0.712820512820513 0.806858974358974;0.00523296703296702 0.726161538461538 0.777869230769231;0.0320450549450548 0.738054945054945 0.747595604395605;0.0882311355311351 0.748630036630037 0.716249084249085;0.140771428571428 0.7584 0.684157142857143;0.175094871794872 0.768316483516483 0.650774725274726;0.200634065934066 0.778392307692308 0.615007692307693;0.230334065934066 0.787759340659341 0.575889010989011;0.273823076923076 0.795176923076923 0.533410989010989;0.329097802197802 0.79999010989011 0.488081318681319;0.387034065934065 0.802712087912088 0.438786813186814;0.45014945054945 0.802179120879121 0.387617582417583;0.517591208791208 0.798130769230769 0.337957142857143;0.583036263736263 0.791689377289377 0.287849084249085;0.646790109890109 0.783124175824176 0.240486813186814;0.708134065934065 0.772743956043956 0.200773626373627;0.766197435897435 0.761349084249084 0.168673992673993;0.820314285714285 0.749814285714286 0.153528571428572;0.869604029304029 0.739435164835165 0.162368864468864;0.914531868131868 0.731723076923077 0.187427472527472;0.954718681318681 0.729240659340659 0.223565934065934;0.986843589743589 0.738140659340659 0.238030402930403;0.996676923076923 0.760668131868132 0.223913186813187;0.995341758241758 0.787452747252747 0.204083516483517;0.988386813186813 0.815494505494505 0.187612087912088;0.976094505494506 0.844460439560439 0.173726373626374;0.965094505494506 0.873561172161172 0.160208058608059;0.960289010989011 0.902178021978021 0.147527472527473;0.961834065934066 0.929990109890109 0.131354945054945;0.968142124542124 0.957045421245421 0.109461538461539;0.9769 0.9839 0.0805000000000006],...
    'FontName','sansserif','FontWeight','bold','Layer','top','XTick',...
    [-70.5 -70.3 -70.1],'XTickLabel',{'-70.5�','-70.3�','-70.1�'},'YTick',...
    [42 42.1 42.2 42.3],'YTickLabel',{'42.0�','42.1�','42.2�','42.3�'});
% Create colorbar
colorbar('peer',subplot2,'FontWeight','bold','FontName','sansserif');

