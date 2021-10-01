Zstate = Z_StateN;

poseID = find(Z_State(:,2)==1);
LMID = find(Z_State(:,2)==2);
count =1;
% %%
% figure
% count=1;
% while count<length(Z_State)
%    if Z_State(i,2)==1
%        plot(Z_State(i,1),Z_State(i+1,1),'k.');
%        count=count+3;
%    else
%        plot(Z_State(count,1),Z_State(count+1,1),'ro');
%        count = count+2;
%    end
%     
% end
%%
poseID = find(Z_State(:,2)==1);
pose_Co = Z_State(poseID,1);
plot(pose_Co(1:3:end,1),pose_Co(2:3:end,1),'k.')
hold on
LMID = find(Z_State(:,2)==2);
LM_Co = Z_State(LMID,1);
plot(LM_Co(1:2:end,1),LM_Co(2:2:end,1),'ro')