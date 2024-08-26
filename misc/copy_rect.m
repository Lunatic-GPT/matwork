function copy_rect()

        h1 = findobj(gca,'Tag','imrect');
        if ~isempty(h1)
          delete(h1);
        end
        h2=findobj(gcf,'Tag','imrect');
        rect_api=iptgetapi(h2);
        imrect(gca,rect_api.getPosition());