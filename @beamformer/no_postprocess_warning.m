function  no_postprocess_warning()
%NO_POSTPROCESS_WARNING To awoid code repetition :)

text = ['You did not call a postprocess, but the channeldata had multiple sequences.',...
          ' This will lead to problems if you dont know what you are doing.'...
          ' Please consider using for example beamforming.go(postprocess.coherent_compound)'];
warning(text);
end

