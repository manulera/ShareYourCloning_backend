import glob
import os
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse, FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from .get_router import get_router

from .endpoints.primer_design import router as primer_design_router
from .endpoints.external_import import router as import_router
from .endpoints.other import router as other_router
from .endpoints.annotation import router as annotation_router
from .endpoints.assembly import router as assembly_router
from .endpoints.no_assembly import router as no_assembly_router
from .endpoints.no_input import router as no_input_router


# ENV variables ========================================
SERVE_FRONTEND = os.environ['SERVE_FRONTEND'] == '1' if 'SERVE_FRONTEND' in os.environ else False
BATCH_CLONING = os.environ['BATCH_CLONING'] == '1' if 'BATCH_CLONING' in os.environ else True

origins = []
if os.environ.get('ALLOWED_ORIGINS') is not None:
    origins = os.environ['ALLOWED_ORIGINS'].split(',')
elif not SERVE_FRONTEND:
    # Default to the yarn start frontend url and the cypress testing
    origins = ['http://localhost:3000', 'http://localhost:5173']

# =====================================================

# Instance of the API object
app = FastAPI()

router = get_router()
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
    expose_headers=['x-warning'],
)


# Workaround for internal server errors: https://github.com/tiangolo/fastapi/discussions/7847#discussioncomment-5144709
@app.exception_handler(500)
async def custom_http_exception_handler(request: Request, exc: Exception):

    response = JSONResponse(content={'message': 'internal server error'}, status_code=500)

    origin = request.headers.get('origin')

    if origin:
        # Have the middleware do the heavy lifting for us to parse
        # all the config, then update our response headers
        cors = CORSMiddleware(
            app=app, allow_origins=origins, allow_credentials=True, allow_methods=['*'], allow_headers=['*']
        )

        # Logic directly from Starlette's CORSMiddleware:
        # https://github.com/encode/starlette/blob/master/starlette/middleware/cors.py#L152

        response.headers.update(cors.simple_headers)
        has_cookie = 'cookie' in request.headers

        print('>>', cors.simple_headers)

        # If request includes any cookie headers, then we must respond
        # with the specific origin instead of '*'.
        if cors.allow_all_origins and has_cookie:
            response.headers['Access-Control-Allow-Origin'] = origin

        # If we only allow specific origins, then we have to mirror back
        # the Origin header in the response.
        elif not cors.allow_all_origins and cors.is_allowed_origin(origin=origin):
            response.headers['Access-Control-Allow-Origin'] = origin
            response.headers.add_vary_header('Origin')

    print(response.headers)
    print(response.body)
    print(response)

    return response


if not SERVE_FRONTEND:

    @router.get('/')
    async def greeting(request: Request):
        html_content = """
            <html>
                <head>
                    <title>Welcome to ShareYourCloning API</title>
                </head>
                <body>
                    <h1>Welcome to ShareYourCloning API</h1>
                    <p>You can access the endpoints documentation <a href="./docs">here</a></p>
                </body>
            </html>
            """
        return HTMLResponse(content=html_content, status_code=200)

else:
    app.mount('/assets', StaticFiles(directory='frontend/assets'), name='assets')
    app.mount('/examples', StaticFiles(directory='frontend/examples'), name='examples')

    @router.get('/')
    async def get_frontend_index(request: Request):
        return FileResponse('frontend/index.html')

    frontend_files = (
        glob.glob('frontend/*.json')
        + glob.glob('frontend/*.ico')
        + glob.glob('frontend/*.png')
        + glob.glob('frontend/*.txt')
    )
    frontend_files = [f.split('/')[-1] for f in frontend_files]

    @router.get('/{name:path}')
    async def get_other_frontend_files(name: str):
        """Catch-all for frontend files"""
        if name in frontend_files:
            return FileResponse(f'frontend/{name}')
        raise HTTPException(404)


app.include_router(primer_design_router)
app.include_router(import_router)
app.include_router(other_router)
app.include_router(annotation_router)
app.include_router(assembly_router)
app.include_router(no_assembly_router)
app.include_router(no_input_router)

if BATCH_CLONING:
    from .batch_cloning import router as batch_cloning_router
    from .batch_cloning.ziqiang_et_al2024 import router as ziqiang_et_al2024_router
    from .batch_cloning.pombe import router as pombe_router

    app.include_router(batch_cloning_router)
    app.include_router(ziqiang_et_al2024_router)
    app.include_router(pombe_router)


# This router must be added last because when SERVE_FRONTEND is True,
# it contains a catch-all route. The catch-all route '/{name:path}' matches any URL path, so if this
# section were placed earlier, it would intercept all requests before they could reach their intended
# API endpoints. For example, requests to '/docs' or '/version' would incorrectly return 404 errors
# instead of reaching their proper handlers.
app.include_router(router)
