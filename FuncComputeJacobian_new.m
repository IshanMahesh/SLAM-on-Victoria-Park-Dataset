function [FX0 J0] = FuncComputeJacobian_new(Zstate,Xstate)

Method=1; % 1 new, 0 or other, old

% compute the F(X0) and Jacobian J0
count=1;count2=1;

if Method==1
    Vector_row = [];
    Vector_column = [];
    Vector_value = [];
else    
    J0=sparse(size(Zstate,1),size(Xstate,1));    
end

temptype=0;tempodoNum=0;tempJodo2=0;tempJLM=0;

while (count<=size(Zstate,1))
    temptype=Zstate(count,2);
    if temptype==1;
        tempodoNum=Zstate(count,3) ;
        [tempJodo1,tempJodo2]=XStateSearch(Xstate,tempodoNum,0,2);
        if tempodoNum==1
            for j=0:2
                if Method==1
                    Vector_row = [Vector_row;count+j];
                    Vector_column = [Vector_column;tempJodo2+j];
                    Vector_value = [Vector_value;1];
                else
                    J0(count+j,tempJodo2+j)=1;
                end
            end
            FX0(count:count+2,1)=Xstate(tempJodo2:tempJodo2+2,1);
        else
            J_small_3=DiffFunction(Xstate(tempJodo1:tempJodo1+2,1),Xstate(tempJodo2:tempJodo2+2,1),1,1);
%            pause
            if Method==1
                Vector_row = [Vector_row;count;count;count;count+1;count+1;count+1;count+2;count+2;count+2];
                Vector_column = [Vector_column;tempJodo1;tempJodo1+1;tempJodo1+2;tempJodo1;tempJodo1+1;tempJodo1+2;tempJodo1;tempJodo1+1;tempJodo1+2];
                Vector_value = [Vector_value;J_small_3(1,1);J_small_3(1,2);J_small_3(1,3);J_small_3(2,1);J_small_3(2,2);J_small_3(2,3);J_small_3(3,1);J_small_3(3,2);J_small_3(3,3)];
            else                
                J0(count:count+2,tempJodo1:tempJodo1+2)=J_small_3;
            end
            
            
            J_small_3=DiffFunction(Xstate(tempJodo1:tempJodo1+2,1),Xstate(tempJodo2:tempJodo2+2,1),2,1);
            
%            pause
            
            if Method==1
                Vector_row = [Vector_row;count;count;count;count+1;count+1;count+1;count+2;count+2;count+2];
                Vector_column = [Vector_column;tempJodo2;tempJodo2+1;tempJodo2+2;tempJodo2;tempJodo2+1;tempJodo2+2;tempJodo2;tempJodo2+1;tempJodo2+2];
                Vector_value = [Vector_value;J_small_3(1,1);J_small_3(1,2);J_small_3(1,3);J_small_3(2,1);J_small_3(2,2);J_small_3(2,3);J_small_3(3,1);J_small_3(3,2);J_small_3(3,3)];                                                       
            else
                J0(count:count+2,tempJodo2:tempJodo2+2)=J_small_3;
            end
            
            c=cos(Xstate(tempJodo1+2,1));s=sin(Xstate(tempJodo1+2,1));
            %F3
            
            FX0(count,1) = (-Xstate(tempJodo1,1)*c-s*Xstate(tempJodo1+1,1)+s*Xstate(tempJodo2+1,1)+Xstate(tempJodo2,1)*c);              %F1
            FX0(count+1,1) = (Xstate(tempJodo1,1)*s-c*Xstate(tempJodo1+1,1)+c*Xstate(tempJodo2+1,1)-Xstate(tempJodo2,1)*s);
            FX0(count+2,1)=Xstate(tempJodo2+2,1)-Xstate(tempJodo1+2,1);
            
        end
        count=count+3;
    else
        tempLMnum=Zstate(count,3);
        tempodoNum=Zstate(count,4);
        [tempJodo1,tempJLM]=XStateSearch(Xstate,tempodoNum,tempLMnum,1);
        if tempodoNum==0
            for j=0:1
                
                if Method==1
                    Vector_row = [Vector_row;count+j];
                    Vector_column = [Vector_column;tempJLM+j];
                    Vector_value = [Vector_value;1];
                else                  
                    
                    J0(count+j,tempJLM+j)=1;
                end
            end
            FX0(count:count+2,1)=Xstate(tempJLM:tempJLM+2,1);
        else
            
            J_small_2=DiffFunction(Xstate(tempJodo1:tempJodo1+2,1),Xstate(tempJLM:tempJLM+1,1),1,2);
            
            %   pause
            
            if Method==1
                Vector_row = [Vector_row;count;count;count;count+1;count+1;count+1];
                Vector_column = [Vector_column;tempJodo1;tempJodo1+1;tempJodo1+2;tempJodo1;tempJodo1+1;tempJodo1+2];
                Vector_value = [Vector_value;J_small_2(1,1);J_small_2(1,2);J_small_2(1,3);J_small_2(2,1);J_small_2(2,2);J_small_2(2,3)];                               
            else
                J0(count:count+1,tempJodo1:tempJodo1+2)=J_small_2;
            end
            
            J_small_2=DiffFunction(Xstate(tempJodo1:tempJodo1+2,1),Xstate(tempJLM:tempJLM+1,1),2,2);
            
            %          pause
            if Method==1
                Vector_row = [Vector_row;count;count;count+1;count+1];
                Vector_column = [Vector_column;tempJLM;tempJLM+1;tempJLM;tempJLM+1];
                Vector_value = [Vector_value;J_small_2(1,1);J_small_2(1,2);J_small_2(2,1);J_small_2(2,2)];               
            else
                J0(count:count+1,tempJLM:tempJLM+1)=J_small_2;
            end
            c=cos(Xstate(tempJodo1+2,1));s=sin(Xstate(tempJodo1+2,1));
            FX0(count,1) = (-Xstate(tempJodo1,1)*c-s*Xstate(tempJodo1+1,1)+s*Xstate(tempJLM+1,1)+Xstate(tempJLM,1)*c);              %F1
            FX0(count+1,1) = (Xstate(tempJodo1,1)*s-c*Xstate(tempJodo1+1,1)+c*Xstate(tempJLM+1,1)-Xstate(tempJLM,1)*s);
            
        end
        count=count+2;
    end
   %   sprintf('Creating Jacobian %.0f %%', (count/size(Zstate,1))*100)
end

if Method==1
%     Vector_row
%     Vector_column
%     pause
    J0=sparse(Vector_row,Vector_column,Vector_value,size(Zstate,1),size(Xstate,1));
end

%spy(J0)

% if Method==1
%     J0_new=J0;
%     save J0_new
% else
%     J0_old=J0;
%     save J0_old
% end
%pause
end
