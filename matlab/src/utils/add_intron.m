%warning: initial commented out section contained updating splicegraph{1}
%as well... check if it makes a difference

%%% adds new introns into splicegraph between idx1 and idx2
%%% if flag1, all end terminal exons in idx1 are preserved
%%% if flag2, all start terminal exons in idx2 are preserved

function splicegraph = add_intron(splicegraph, idx1, flag1, idx2, flag2)

if(nargin>3),
    adj_mat = triu(splicegraph{2}) ;

    if flag1>0,
        for i1=idx1,
            %%% if exon is end-terminal
            if all(adj_mat(i1,:)==0),
                splicegraph{1}(:,end+1) = splicegraph{1}(:, i1) ;
                splicegraph{2}(end+1, end+1) = 0 ;
                splicegraph{2}(:,end) = splicegraph{2}(:, i1) ;
                splicegraph{2}(end,:) = splicegraph{2}(i1, :) ;
				if length(splicegraph)>2
					splicegraph{3}(:, end+1) = splicegraph{3}(:, i1);
				end
            end ;
        end ;
    end ;

    if flag2>0,
        for i2=idx2,
            %%% if exon is start-terminal
            if all(adj_mat(:,i2)==0),
                splicegraph{1}(:,end+1) = splicegraph{1}(:, i2) ;
                splicegraph{2}(end+1, end+1) = 0 ;
                splicegraph{2}(:,end) = splicegraph{2}(:, i2) ;
                splicegraph{2}(end,:) = splicegraph{2}(i2, :) ;
				if length(splicegraph)>2
					splicegraph{3}(:, end+1) = splicegraph{3}(:, i2);
				end
            end ;
        end ;
    end ;
end

for i1=idx1,
    for i2=idx2,
        splicegraph{2}(i1,i2) = 1 ;
        splicegraph{2}(i2,i1) = 1 ;
    end ;
end ;

