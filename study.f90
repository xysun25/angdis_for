module atomdef   ! ��װatomdef����
    implicit none
    type atom    ! �Զ���һ����atom����������,�����������������������������Ϊatom
        character(len=10)::element   ! ��¼Ԫ������
        real::x,y,z,m,e,sig,eps   ! ��¼ԭ�������� L-J ���нضϾ��������
        end type   ! �Զ����������ͽ���
    
    type mol   ! �Զ���һ����mol ����������
        character(len=10)::numberofmol   ! ��¼������Ŀ
        real::x,y,z   ! ��¼��������
        end type   ! �Զ����������ͽ���
    
end module 

! ������
program main
    use atomdef   ! ����atomdefģ��
    implicit none
    character(len=20)::filename,mmol1,mmol2,mmol3,trjfile   ! ��¼�ļ�����
    integer::num, numberofAIM
    integer i,j,k,u,l,p,q,r,trj,ls,u1,u2,u3,ff1,ff2,lnn,np,nn,mm,up1,upan,uu1,uu2,uu3,nsite,tt,qt
    integer numberoftrj,o1,h1,oo,ohh,ha
    integer num1,num2,num3,numberofAIM1,numberofAIM2,numberofAIM3
    integer numberofatom1,numberofatom2,numberofatom3,numberofatoms_2
    type(atom),allocatable::atoms_1(:),atoms_2(:),atoms(:,:)   ! ��������atom���͵ı������������ɱ�һά���飬��һ��atom���͵ı������ɱ��ά����
    real X,Y,Z,M,Xcom,Ycom,Zcom,x00,y00,dr,high,oh,ohx,ohy,ohz,ohl,oohx,oohy,oohz,oohhx,oohhy,oohhz,oohhl,ohhx,ohhy,ohhz,ool
    real a0,b0,c0,a,b,c,ax,by,cz,interv,theta,bin,angle,c1,c2,thetah,an,bn,cn,anx,any,anz,bnx,bny,bnz,cnx,cny,cnz
    integer numberofinterv
    integer::numberofions,numberofatoms
    real::o(2),h(2)
    real,allocatable::Hang(:),angdis(:,:,:)   ! ����һ���ɱ�һά���飬һ���ɱ���ά����
    
    ! Ҫдin�ļ���Ϊ���룬����������������Ͳ����ļ�������
    read (*,*) filename   ! in�ļ������������ļ����ƺͱ�������
    write(*,*) 'the input file is:', filename
    open(100,file=filename,status="old")
    read(100,*) trjfile,numberoftrj
    write(*,*) 'trjfile and numberof trj',trjfile,numberoftrj
    read (100,*) mmol1,mmol2,mmol3
    write(*,*) 'the mmol files are:', mmol1,mmol2,mmol3
    read (100,*) num1,num2,num3
    write(*,*) 'the number of molecules are:', num1,num2,num3
    read(100,*) a0,b0,c0,ax,by,cz
    read(100,*) interv, numberofinterv, bin
    read(100,*) u1, u2, u3
    read(100,*) uu1, uu2, uu3
    
    o=[4,12]
    h=[16,18]
    o1=19
    h1=20
    
    ! ����EmimOH.mmol �ļ�
    open(11,file=mmol1,status="old")   ! ���Ѵ��ڵ�mmol1�ļ���ָ���ļ�����Ϊ11
    read(11,*)   ! ������ļ�
    read(11,*) numberofAIM1    ! 20
    allocate(atoms_1(numberofAIM1))   ! allocate ����atoms_1һά������ڴ�ռ䣬����Ҫ�ͷŵ�
    ! ѭ��
    do i=1,numberofAIM1   ! i��ʼ��ֵ��numberofAIM1��ֹ��ֵ������Ĭ��Ϊ1;ѭ��20��
        ! ʹ��atom���ͱ���atoms_1��element��Ԫ�صķ�����������Ԫ����%���
read(11,*)atoms_1(i)%element,atoms_1(i)%x,atoms_1(i)%y,atoms_1(i)%z,atoms_1(i)%m,atoms_1(i)%e,atoms_1(i)%sig,atoms_1(i)%eps        
    end do
    rewind(11)
    
    ! ����NTf2.mmol �ļ�
    open(12,file=mmol2,status="old")
    read(12,*)
    read(12,*) numberofAIM2   ! 14
    allocate (atoms_2(numberofAIM2))
    do j=1,numberofAIM2
read(12,*)atoms_2(j)%element,atoms_2(j)%x,atoms_2(j)%y,atoms_2(j)%z,atoms_2(j)%m,atoms_2(j)%e,atoms_2(j)%sig,atoms_2(j)%eps
    end do
    rewind(12)
    
    ! ����mxene.mmol �ļ�
    open(13,file=mmol3,status="old")
    read(13,*) numberofAIM3   ! 14
    close(13)
    
    numberofatom3=num3*numberofAIM3   ! 600*14=8400
    num=num1*numberofAIM1+num2*numberofAIM2+num3*numberofAIM3   ! 500*20+500*15+600*14=25900
    numberofAIM=numberofAIM1+numberofAIM2   ! 20+15=35
    nsite=int(180/bin)   ! 180/0.5=360
    mm=num1*num3   ! 500*600=300000
    ha=2/numberofinterv   !2/100(140)(70)
     
    ! �򿪹켣�ļ�
    open(10,file=trjfile,status="old")
    allocate(angdis(4,numberofinterv,nsite))   ! ������ά�ɱ�����angdis���ڴ�angdis(4,100,360)
    allocate(Hang(nsite))   ! ����һά����Hang���ڴ�Hang(360)
    
    open(16,file="cation_angl.dat",status="unknown")
    open(17,file="cation_ang2.dat",status="unknown")
    open(18,file="anion_angl.dat",status="unknown")
    open(19,file="anion_ang2.dat",status="unknown")
    open(20,file="hydrogen_ang.dat",status="unknown")
    
    ff1=0
    ff2=0
    upan=0
    up1=0
    np=100
    nn=0
    ! �ڹ켣�ļ���ѭ��
    do trj=1,numberoftrj   ! trj��ʼ��ֵ��numberoftrj��ֹ��ֵ��֡������
        allocate(atoms(num,numberofAIM))   ! atoms���Զ�����������atom��һ����ά����atoms(25900,35)
        read(10,*)
        read(10,*)
        do k=1,num1   ! ѭ��500��
            do u=1,numberofAIM1   ! ѭ��20��
                read(10,*) atoms(k,u)%element,atoms(k,u)%x,atoms(k,u)%y,atoms(k,u)%z
            end do
        end do 
        do k=num1+1,num2+num1
		    do u=1,numberofAIM2
			    read (10,*)atoms(k,u)%element,atoms(k,u)%x,atoms(k,u)%y,atoms(k,u)%z
		    end do
	    end do

	    do k=num2+num1+1,num2+num1+num3
		    do r=1,numberofAIM3
			    read(10,*) atoms(k,r)%element,atoms(k,r)%x,atoms(k,r)%y,atoms(k,r)%z
		    end do
        end do


    

        
    
            
    
    
    
    
    
    