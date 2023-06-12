function [ majorVersion,minorVersion,patchVersion ] = extractVersionInfo( versionString )
% Extracts major, minor and patch version number from any version string.
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) embotech AG, Zurich, Switzerland, 2013-2023.

    versionToken = regexp(versionString,'\w*(\d+\.\d+(?:\.\d+)?)\w*','tokens');
    if (length(versionToken) >= 1)
        versionInfo = strsplit( versionToken{1}{1},'.' );
        if (length(versionInfo) >= 1)
            majorVersion = str2double(versionInfo{1});
        else
            majorVersion = NaN;
        end
        if (length(versionInfo) >= 2)
            minorVersion = str2double(versionInfo{2});
        else
            minorVersion = NaN;
        end
        if (length(versionInfo) >= 3)
            patchVersion = str2double(versionInfo{3});
        else
            patchVersion = NaN;
        end
    else
        majorVersion = NaN;
        minorVersion = NaN;
        patchVersion = NaN;
    end
end
