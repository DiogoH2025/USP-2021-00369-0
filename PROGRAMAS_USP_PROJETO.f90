!############################################################
!Criação de Redes obedecendo o algoritmo UCM
subroutine network
implicit none

!=================================
!grau de corte UCM
	kc=2.d0*((1.d0*N)**(0.5d0))
!=================================
	allocate(pk(kc))
!=================================
!norma da distribuicao e de knn
	norma=0.d0
	do grau=ko,kc

		pk(grau)=(1.d0*grau)**(-1.d0*gama)
		norma=norma+(1.d0*grau)**(-1.d0*gama)

	enddo
	
	norma=1.d0/norma
	pk=norma*pk
!=================================
!	distribuibuição de grau
1	xo=(1.d0*ko-1.d0)**(1.d0-gama)
	xc=(1.d0*kc)**(1.d0-gama)
	expo=1.d0/(1.d0-gama)

	nedge=0
	do cont=1,N
		do
			iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
			x=(xo-z*(xo-xc))**expo
			grau=ceiling(x)
			iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
			if((z*norma/x**gama) >= pk(grau))cycle
			nedge=nedge+grau
			degree(cont)=grau
			exit
		enddo
	enddo

	if(mod(nedge,2)/=0) go to 1

	allocate(lista(nedge),matriz(nedge))
!=================================
2	begin=0;nedge=0
	kmed=0.d0;k2med=0.d0
	do cont=1,N
		begin(cont)=nedge
		do ind1=1,degree(cont)
			nedge=nedge+1;lista(nedge)=cont
		enddo
		kmed=kmed+1.d0*degree(cont)
		k2med=k2med+1.d0*degree(cont)*degree(cont)
	enddo
!=================================
	aux=0;matriz=0
	do

		if(nedge==0)exit

		ntent=0
		do

			iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
			ind1=z*nedge+1!posição na lista
	!=================================
			iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
			ind2=z*nedge+1!posição na lista
	!=================================
			if(ntent==try_lim)go to 2
	!=================================
			ntent=ntent+1
			if(lista(ind1)==lista(ind2))cycle !evita autoconexões
			flag=0		
			do cont=1,aux(lista(ind1)) 
				if(matriz(begin(lista(ind1))+cont)==lista(ind2))then
					flag=1;exit
				endif
			enddo
			if(flag==1)cycle !evita múltiplas conexões
			exit
		enddo
	!=================================
	!vértice lista(i)
		aux(lista(ind1))=aux(lista(ind1))+1
		matriz(begin(lista(ind1))+aux(lista(ind1)))=lista(ind2)			
	!vértice lista(j)
		aux(lista(ind2))=aux(lista(ind2))+1
		matriz(begin(lista(ind2))+aux(lista(ind2)))=lista(ind1)			
	!ajuste lista
		if(ind1<ind2)then
			cont=ind1;ind1=ind2;ind2=cont
		endif
		lista(ind1)=lista(nedge);nedge=nedge-1
		lista(ind2)=lista(nedge);nedge=nedge-1
	enddo

	deallocate(pk)
!=================================

end subroutine network
!#########################################################################
!Implementação do algoritmo de Gillespie para Dinâmica SIRS
subroutine SIRS_multiple_stages
implicit none
	
	flag=0
	tempo=0.d0;nconf=0
	Nciclo=0;dTIME=0	
	norma=0.d0;pqst=0.d0
	pqst_susc=0.d0;pqst_rec=0.d0
	phi=0.d0
		
	call condicao_inicial

	do
		
		!escolho o tipo de processo   	
		processo(1)=(mu*nstage*Ninf)/(lambda*Ne+mu*nstage*Ninf+alpha*Nrec)	! evolução do quadro/cura		
		processo(2)=(alpha*Nrec)/(lambda*Ne+mu*nstage*Ninf+alpha*Nrec)	! returno ao suscetível 		
		processo(3)=(lambda*Ne)/(lambda*Ne+mu*nstage*Ninf+alpha*Nrec)		! infecção

		iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
		prob=0.d0			
		do nproc=1,3
			prob=prob+processo(nproc)
			if(z<prob)exit
		enddo
		
		select case(nproc)
		!====================================
		!cura
			case (1) 
			
				!------------------------------
				!lista de disposição do número de infectados
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
				ind=Ninf*z+1; vert=inf(ind)	
				stage=rede(vert)				
				if(stage/=nstage)then  !evolução do quadro
					rede(vert)=stage+1
				else !	processo de cura		
											
					rede(vert)=-1;Ne=Ne-degree(vert)				
					inf(ind)=inf(Ninf);Ninf=Ninf-1		
					Nrec=Nrec+1;recover(Nrec)=vert
				!------------------------------
					interval=tempo-t0_inf(vert)!;t0_inf(vert)=0.d0
					ncaixa=int(interval/delta)
					dTIME(ncaixa)=dTIME(ncaixa)+1		
					Nciclo=Nciclo+1						
				!------------------------------	
					if(Ninf==0)call quasi_estacionario
					
				endif
		!===================================
		!retorno para suscetível	
			case(2)
				
				!lista de recuperados
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0				
				ind=z*Nrec+1;vert=recover(ind)
				recover(ind)=recover(Nrec);Nrec=Nrec-1
				rede(vert)=0
		!===================================			
		!infecção de suscetível
			case default
			
				!sorteio de um infectado proporcionalmente ao grau
				do
					iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
					ind=nvert*z+1
					vert=lista(ind)
					if((rede(vert)/=0).and.(rede(vert)/=-1)) exit		
					aux(vert)=aux(vert)-1
					lista(ind)=lista(nvert);nvert=nvert-1
				enddo
				!-------------------------------
				!escolha de um vizinho aleatoriamente
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
				iviz=z*degree(vert)+1
				viz=matriz(begin(vert)+iviz)					
				if(rede(viz)==0)then
					rede(viz)=1;Ne=Ne+degree(viz)
					Ninf=Ninf+1;inf(Ninf)=viz
				!------------------------------					
					t0_inf(viz)=tempo !distribuição de tempo de infecção
				!------------------------------
				!ajusta a lista
					dif=degree(viz)-aux(viz)
					aux(viz)=degree(viz)
					do cont=1,dif
						nvert=nvert+1
						lista(nvert)=viz
					enddo
				!-------------------------------
				endif
			
		end select 
	!========================================	
		dt=1.d0/(lambda*Ne+mu*nstage*Ninf+alpha*Nrec) 
		tempo=tempo+dt
		call historico
		
		if(tempo>=trlx)then
		
			pqst(Ninf)=pqst(Ninf)+dt
			pqst_rec(Nrec)=pqst_rec(Nrec)+dt
			pqst_susc(N-Nrec-Ninf)=pqst_susc(N-Nrec-Ninf)+dt
			norma=norma+dt
			
		!----vetor de localização-----	
			do ind=1,Ninf
				vert=inf(ind)
				phi(vert)=phi(vert)+dt
			end do
	
			if(flag==0)then
				t0=tempo;flag=1
			endif
			
		endif
!======================================
!reesquentar o número aleatório aqui 
!======================================
		if(tempo>=tf)then
		!------------------------------
			do vert=1,N		
				if( (rede(vert)==0).or.(rede(vert)==-1) )cycle
				interval=tempo-t0_inf(vert);t0_inf(vert)=0.d0
			enddo										
			exit
		endif
!======================================		
		
	enddo

	pqst=pqst/norma
	pqst_rec=pqst_rec/norma
	pqst_susc=pqst_susc/norma	
	call estatistica
!-------------------------------		
	phi=phi/(tempo-t0)	
	call distribuicao_phi	
	
end subroutine SIRS_multiple_stages
!#########################################################################
!Estabelece a condição inicial da dinâmica (todos infectados) 
subroutine condicao_inicial
implicit none

	integer::vert,grau

	Ninf=0;Ne=0;Nrec=0
	rede=0;aux=0;nvert=0
	
	do vert=1,N
		
		rede(vert)=1;t0_inf(vert)=tempo
		Ninf=Ninf+1;inf(vert)=vert
		Ne=Ne+degree(vert);aux(vert)=degree(vert)
			
		do grau=1,degree(vert)
			nvert=nvert+1;lista(nvert)=vert
		enddo

	enddo
	
end subroutine condicao_inicial
!#########################################################################
!Rotina que recupera uma das configurações passadas na dinâmica
subroutine quasi_estacionario
implicit none

	integer::ind,cont,conf,vert

	iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
	conf=z*nconf+1	
	Ninf=sick(conf);Nrec=treat(conf)
!----------infectado---------------
	do ind=1,nvert
		aux(lista(ind))=0
	enddo

	nvert=0;Ne=0;rede=0
	do ind=1,Ninf
	
		inf(ind)=hist(conf,ind)
		vert=inf(ind)
		rede(vert)=Istage(conf,ind)
		t0_inf(vert)=tempo
		
		do cont=1,degree(vert)
			nvert=nvert+1
			lista(nvert)=vert
		enddo
		
		aux(vert)=degree(vert)
		Ne=Ne+degree(vert)
		
	enddo
!---------recuperado----------------
	do ind=1,Nrec
		recover(ind)=hist(conf,Ninf+ind)
		vert=recover(ind)
		rede(vert)=-1		
	enddo
	
end subroutine quasi_estacionario
!#########################################################################
!Rotina que permite adicionar tubos com probabilidade 
!proporcional ao grau de um veŕtice núcleo.
subroutine LINK_TUBE_SUPER_PREF
	implicit none

!============================================
!Esquenta o gerador de número aleatório
	do vert2=N+1,N+N1
	
	!---------------------------
		aux(vert2)=0;kgrau(vert2)=0
	!---------------------------
		nviz=0
		do
					
			do
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
				vert1=z*N+1!posição na lista
				plink=((1.d0*degree(vert1))/(1.d0*KHUB))**2.d0
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0				
				if(z<plink)exit
			enddo	 			
		!---------------------------
		!evita múltiplas conexões
			flag=0
			do ind1=1,nviz
				if(viz(ind1)==vert1) flag=1 
			enddo
			if(flag==1)cycle
		!--------------------------		
			aux(vert1)=aux(vert1)+1
			aux(vert2)=aux(vert2)+1
			nviz=nviz+1;viz(nviz)=vert1
			if(nviz==2)exit
		!--------------------------		
		enddo
		
	enddo
!============================================
	nlink=0	
	do vert1=1,N+N1
		
		start(vert1)=nlink
		nlink=nlink+aux(vert1)
				
	enddo 
!-----------------------------	
	allocate(Adj(nlink))
!-----------------------------
	do vert1=1,N
	
		kgrau(vert1)=degree(vert1)		
		do grau=1,degree(vert1)
			Adj(start(vert1)+grau)=matriz(begin(vert1)+grau)
		enddo
		
	enddo 	
!===========================================
!repetição/construção matrix
	call esquenta			
	do vert2=N+1,N+N1
		
		nviz=0
		do
							
			do
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
				vert1=z*N+1!posição na lista
				plink=((1.d0*degree(vert1))/(1.d0*KHUB))**2.d0
				iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0				
				if(z<plink)exit
			enddo	 			
		!---------------------------
		!evita múltiplas conexões
			flag=0
			do ind1=1,nviz
				if(viz(ind1)==vert1) flag=1 
			enddo
			if(flag==1)cycle
		!--------------------------			
			kgrau(vert1)=kgrau(vert1)+1
			kgrau(vert2)=kgrau(vert2)+1
		!---------------------------		
		!vértice lista(i)
			Adj(start(vert1)+kgrau(vert1))=vert2			
		!vértice lista(j)		
			Adj(start(vert2)+kgrau(vert2))=vert1		
		!---------------------------					
			nviz=nviz+1;viz(nviz)=vert1			
			if(nviz==2)exit
		!--------------------------		
		enddo		
				
	enddo
!===========================================
	
end subroutine LINK_TUBE_SUPER_PREF
!###########################################################
!Rotina para criar lista de configurações a serem recuperados
!quando o sistema atinge o estado absorvente
subroutine historico
implicit none

	if(nconf==ntconf)then

		iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
		conf=z*ntconf+1
		
		iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
		if(z<pre*dt)then
		!-------infectados-----------
			sick(conf)=Ninf
			do ind=1,Ninf
				vert=inf(ind)
				hist(conf,ind)=vert
				Istage(conf,ind)=rede(vert)
			enddo		
		!------recuperados-----------	
			treat(conf)=Nrec
			do ind=1,Nrec
				hist(conf,Ninf+ind)=recover(ind)
			enddo
		!----------------------------									
		endif

	else 

		iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
		if(z<dt)then
			nconf=nconf+1
		!------infectados------------
			sick(nconf)=ninf
			do ind=1,ninf
				vert=inf(ind)
				hist(nconf,ind)=vert
				Istage(nconf,ind)=rede(vert)
			enddo
		!------recuperados-----------	
			treat(nconf)=nrec
			do ind=1,nrec
				hist(nconf,ninf+ind)=recover(ind)
			enddo
		!----------------------------		
		endif

	endif

end subroutine historico
!#########################################################################
!Computar o IPR da atividade através da distribuição de atividade
subroutine distribuicao_phi
implicit none

!------------------------------------	
	norma=0.d0
	do vert=1,N
		norma=norma+phi(vert)*phi(vert)			
	enddo
	norma=sqrt(norma)	
	phi=phi/norma	
!-----------------------------------
	ipr=0.d0
	do vert=1,N	
	!================================
		ipr=ipr+phi(vert)**4.d0
	!================================		
	enddo	
!-----------------------------------
		
end subroutine distribuicao_phi
!#########################################################################
!Computa as quantidades quasi-estacionárias
!Suscetibilidade, densidade e lifespan
subroutine estatistica
implicit none

	rho=0.d0;rho2=0.d0
	drec=0.d0;dsus=0.d0
	do ind=1,N
	
		razao=(1.d0*ind)/(1.d0*N)
	!-----infectados------------------
		rho=rho+razao*pqst(ind)
		rho2=rho2+razao*razao*pqst(ind)
	!-----recuperados-----------------
		drec=drec+razao*pqst_rec(ind)
	!-----suscetiveis-----------------
		dsus=dsus+razao*pqst_susc(ind)
		
	enddo

	lifespan=1.d0/pqst(1)
	susc=1.d0*N*((rho2-rho*rho)/rho)

end subroutine estatistica
!#########################################################################
!Diagonalização matriz Jacobiana do SIRS
subroutine metodo_potencia
implicit none

	do ind1=1,N	
		iseed=(AMAG*iseed)+BMAG; rcont=rcont+1; z=r_1*iseed+0.5d0
		vec(ind1)=z
	enddo
	
  eigmax_old=0.d0;eigmax=0.d0
  niter=0;dif=1.d0

  norma1=sum(vec**2.d0)**0.5d0
  eigmax_old=0.d0
  IPR_old=0.d0
  factor=2 !	
	
	fator1=const*lambda*lambda  
	fator2=const*lambda*(mu+lambda+a)
	
	do while(dif>tol.or.niter<=1000) 

		niter=niter+1
		vec=vec/norma1	
		vecaux2=vec

		do mult=1,factor
			vecaux1=vec
			do ind1=1,N
				nviz=degree(ind1);inicio=begin(ind1)
				soma=0.d0
				do ind2=1,nviz
					vert=matriz(inicio+ind2)
					soma=soma+vecaux1(vert)
				enddo	
				vec(ind1)=fator1*(kmax-nviz)*vecaux1(ind1)+fator2*soma
			enddo
		enddo


		norma1=sum(vec**2.d0)**0.5d0
		norma2=sum(vecaux1**2.d0)**0.5d0
		IPR=sum(((vec/norma1+vecaux1/norma2)*0.5d0)**4.d0) !to circunvent the non-monotonic covnergence
		eigmax=dot_product(vec,vecaux2)**(1.d0/factor)
		dif=max(abs(IPR-IPR_old),abs(eigmax-eigmax_old))
		eigmax_old=eigmax
		IPR_old=IPR

	enddo


end subroutine metodo_potencia
!#########################################################################
!Ajustes dos valores para convergência da Diagonalização de Pares
subroutine convergencia_DIAGONALIZACAO_PARES

	tol=1e-7
	lambda=0.01d0
	dlambda=0.005d0

	do

	!valores iniciais 
	!========================
		const=1.d0/(a+mu+2.d0*lambda)
		call metodo_potencia
		jacob=-1.d0*(mu+const*kmax*lambda*lambda)+eigmax		

		if(jacob<0.d0)then
			lambda=lambda+dlambda

		else
			if(jacob<tol)then
				exit
			else
				lambda=lambda-dlambda
				dlambda=dlambda/2.d0
				lambda=lambda+dlambda
			endif
		endif
!========================
	!reesquentar o número aleatório  
!========================
	enddo
	
end subroutine convergencia_DIAGONALIZACAO_PARES
!#########################################################################
!Computar o fator necessário para calcular o MMCA
subroutine fator_qi
	
	do vert1=1,N !vértice $i$
	
		fator=1.d0
		nviz=degree(vert1)

		do ind1=1,nviz
		
			vert2=matrix(begin(vert1)+ind1)		!vértice $j$						
			fator=fator*(1.d0-beta*rho(vert2))
			!fator=fator*(1.d0-beta*(1.d0-(1.d0-wij(vert2)/wi(vert1))**lambda_i)*rho(vert2))!termo geral			
		enddo
		
		qi(vert1)=fator	
		
	enddo
		
end subroutine fator_qi
!=======================================================================
!Computar a densidade de infectados do MMCA
Subroutine evolution_rho

	mrho=0.d0
	do vert1=1,N  ! vértice $i$

		rho(vert1)=(1.d0-qi(vert1))*(1.d0-rho_past(vert1))+(1.d0-mu)*rho_past(vert1)+mu*(1.d0-qi(vert1))*rho_past(vert1)
		mrho=mrho+rho(vert1)
			
	enddo

End Subroutine evolution_rho
!#########################################################################
