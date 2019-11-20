      program scattering
      use input
      use channels 
      use scattwf
      implicit none 
      
      ! reading in the parameters 
      call initialize()
      
      ! check the input and make the local copy 
      call check()
      
      ! setup the channel index 
      call alpha_2b()
      
      ! solve the scattering problem 
      call scatt2b()
      
      ! list the output informations
      call fkind()
      
      
      end program 