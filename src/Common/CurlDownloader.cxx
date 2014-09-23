// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/progress.hpp>

#include <curl/curl.h>

#include "Common/PE.hh"
#include "Common/CurlDownloader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

CurlDownloader::CurlDownloader()
{
}

//////////////////////////////////////////////////////////////////////////////

CurlDownloader::~CurlDownloader()
{
}

//////////////////////////////////////////////////////////////////////////////

inline
int CurlDownloader::write_data(void * buffer, int size, int nmemb, void * data)
{
   class FileHandle *fh = static_cast<FileHandle*>(data);
   if (!fh->isopen)
   {
     fh->fp = fopen(fh->fpath.string().c_str(),"w");
     fh->isopen = true;
   }
   return fwrite(buffer, size, nmemb, fh->fp);
}

//////////////////////////////////////////////////////////////////////////////

int CurlDownloader::progress_func( class FileHandle * fh,
                  double dltotal,
                  double dlnow,
                  double ultotal,
                  double ulnow)
{
   unsigned int now = static_cast< unsigned int >(dlnow);
   if (fh->progress == NULL)
   {
      fh->tot = static_cast< unsigned int >(dltotal);
      fh->progress = new boost::progress_display(fh->tot);
   }

   for (unsigned int i = fh->curr; i < now && i < fh->tot; ++i)
   {
      ++(*(fh->progress));
   }

   fh->curr = now;
   return 0; /* non zero would abort download */
}

//////////////////////////////////////////////////////////////////////////////

void CurlDownloader::download ( const std::string& url, const std::string& filepath )
{
  // downloader works serially on the processor with rank == 0
  if (PE::GetPE().IsParallel())
  {
    // all processes come here
    PE::GetPE().setBarrier();

    // only rank 0 downloads
    if (PE::GetPE().GetRank() == 0)
    {
      do_download(url,filepath);
    }

    // all processes come here
    PE::GetPE().setBarrier();
  }
  else
  {
    do_download(url,filepath);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CurlDownloader::do_download ( const std::string& url, const std::string& filepath )
{
   CURL *curl;
   CURLcode res;

   class FileHandle fhandle;
   fhandle.fpath = filepath;

   curl = curl_easy_init();

   if (curl) {

         // skip SSL verifications
         curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0);
         curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0);
         curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1);

         /* set URL to get */
         curl_easy_setopt(curl, CURLOPT_URL, url.c_str());

         /* set function to write the data */
         curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, CurlDownloader::write_data);

         /* set user local data to pass to write function */
         curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void*)&fhandle);

         curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0);
         curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, CurlDownloader::progress_func);
         curl_easy_setopt(curl, CURLOPT_PROGRESSDATA, &fhandle);

         /* do it */
         res = curl_easy_perform(curl);
         if (res != CURLE_OK)
         {
           std::string msg ( curl_easy_strerror(res) );
           curl_easy_cleanup(curl);
           throw URLException (FromHere(),msg);
         }

         /* always cleanup */
         curl_easy_cleanup(curl);
   }

   if (fhandle.isopen)
   {
     fclose(fhandle.fp);
     fhandle.isopen = false;
   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Common

} // namespace COOLFluiD
