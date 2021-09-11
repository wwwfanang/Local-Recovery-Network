
library(igraph)

N<- 1000
k<- 10  
M<- 4
tau<- 100
tau1<- 1

sep<- 0.02
p_star<- c(seq(0,0.99,sep),seq(0.99,0,-sep))
p<- -log(1-p_star,exp(1))/tau
m<- 4
r<- 0.78

iter<- ITER

##  Net
g<- sample_k_regular(N, k)
V(g)$name<- as.character(1:N)

state<- data.frame(internal= rep(0,N),external= rep(0,N))

z_stable<- vector()
for(i in seq_along(p_star)){
  
  num_z<- 1
  z<- 0
  repeat{
    
    ## Internal Failure - LA
    id<- as.vector(components(g)$membership)==which.max(components(g)$csize)
    v_sam<- sample(V(g)$name[id],1)
    
    num<- 1
    repeat{
      if((ego_size(g,num,v_sam)>=((p[i])*N))||
         ((ego_size(g,num,v_sam)==max(components(g)$csize))))
        break
      num<- num+1
    }
    if(ego_size(g,num,v_sam)>=((p[i])*N)){
      v_att<- unlist(ego(g,num,v_sam))[1:((p[i])*N)]
    }else{
      v_att<- unlist(ego(g,num,v_sam))
      v_att<- c(v_att,sample(V(g)$name[-id],((p[i])*N-length(v_att))))
    }
    state$internal[which(V(g)$name%in%v_att)]<- tau
    
    ## External Failure
    id_f<- union(which(state$internal!=0),which(state$external!=0))
    nei<- ego(g,1,V(g)$name)
    nei_active<- lapply(nei, function(x) setdiff(x$name,unique(c(x$name[1],V(g)$name[id_f]))))
    
    id_nei<- which(lengths(nei_active)<=m)
    id_nei<- setdiff(id_nei,which(state$internal!=0))
    rr<- runif(length(id_nei),0,1)
    id_r<- which(rr<= r)
    state$external[id_nei[id_r]]<- tau1
    
    ## Recovery
    id_tau<- setdiff(which(state$internal>0),which(V(g)$name%in%v_att))
    state$internal[id_tau]<- state$internal[id_tau]-1
    
    id_tau1<- setdiff(which(state$external>0),id_nei[id_r])
    state$external[id_tau1]<- state$external[id_tau1]-1
    
    num_z<- num_z+1
    z[num_z]<- length(intersect(which(state$internal==0),which(state$external==0)))/N
    print(c(p_star[i],num_z,length(which(state$internal==0))/N,length(which(state$external==0))/N,length(intersect(which(state$internal==0),which(state$external==0)))/N))
    if((num_z>=200)&(mean(tail(diff(z)))<0.001))
      break
  }
  z_stable[i]<- tail(z,1)
}
write.table(data.frame(p_star,z_stable),paste('RR-Recover_N=',N,'_k=',k,'_m=',m,'_tau=',tau,'_tau1=',tau1,'_r=',r,'_iter=',iter,'_t-z.txt',sep=''),row.names = F)
