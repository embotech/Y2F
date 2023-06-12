function [ answer ] = satisfiesMinimalVersion( responseString,minVersionString )
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) embotech AG, Zurich, Switzerland, 2013-2023.

    if (nargin < 2)
        minVersionString = '0.0.0';
    end

    [ majorVersion,minorVersion,patchVersion ] = extractVersionInfo( responseString );
    [ minMajorVersion,minMinorVersion,minPatchVersion ] = extractVersionInfo( minVersionString );

    if (majorVersion > minMajorVersion) || ...  
       ((majorVersion == minMajorVersion) && (minorVersion > minMinorVersion)) || ...
       ((majorVersion == minMajorVersion) && (minorVersion == minMinorVersion) && (patchVersion >= minPatchVersion))
        answer = true;
    else
        answer = false;
    end
end

