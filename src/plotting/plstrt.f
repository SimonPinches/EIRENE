        subroutine eirene_plstrt
        external grstrt
cdr  initialize proprietary EIRENE "GR plotting software"
        call grstrt(35,8)
        return
        end subroutine eirene_plstrt
