
1) Deploy Filemage: https://www.filemage.io/docs/gcp.html
2) Create Service Account, set permissions as defined above (ensure IAM Service Account Credentials API is enabled)
  - create SA, and add a custom role with both `Storage Object Admin` & `Service Account Token Creator` permisson sets
3) create GCP bucket as endpoint for FileMage
4) login & update password: https://www.filemage.io/docs/quickstart.html
5) create key for service account to enable endpoint access
5) create endpoint: https://www.filemage.io/docs/endpoints.html
6) create users: https://www.filemage.io/docs/users.html
7) enable workspace portal
- add the following to `/etc/filemage/config.yml` (must use `sudo`)
```
management_port: 8443
workspace_port: 443
```
8) enable CORS on all storage buckets: https://www.filemage.io/docs/cors.html#configure-cors-for-azure-blob-storage
- Create the following file, specific for your use case
```
[
  {
    "maxAgeSeconds": 3600,
    "method": ["PUT"],
    "origin": [
      "https://34.57.66.164"
    ],
    "responseHeader": ["content-type"]
  }
]
```
- `gcloud storage buckets update gs://<BUCKET> --cors-file=<CORS_FILE>` to update cors setup
- `gcloud storage buckets describe gs://<BUCKET>` to describe current cors setup (`cors_config` block)
