      function tunedat(x,tuned,at,iopt,number)
c   iopt=0   interpolation in linear scale
c   iopt=1   interpolation in logarithmic scale
c Modified by Maximilien Fechner, march 2007
c number is 5 (Sk-1 and 2) or 4 (SK-III)

      integer number
      real tuned(5),at(5)

      if(x.lt.at(1))then
        tunedat=tuned(1)
      else
        do i=2,number
          if(x.lt.at(i))then
            if(iopt.eq.0)then
              tunedat=(tuned(i-1)*(at(i)-x)+tuned(i)*(x-at(i-1)))
     &                   /(at(i)-at(i-1))
            else
              tunedat=exp((log(tuned(i-1))*(at(i)-x)
     &                    +log(tuned(i))*(x-at(i-1)))
     &                   /(at(i)-at(i-1)))
            endif
            goto 10
          endif
        enddo
        tunedat=tuned(number)
      endif
   10 continue
      return
      end
