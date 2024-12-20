import unittest
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.testclient import TestClient

import shareyourcloning.main as _main
from shareyourcloning.api_config_utils import custom_http_exception_handler


client = TestClient(_main.app)


class GreetingTest(unittest.TestCase):
    def test_greeting(self):
        response = client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<a href="./docs">here</a>', response.text)


class InternalServerErrorTest(unittest.IsolatedAsyncioTestCase):

    async def test_internal_server_error(self):
        request = Request(scope={'type': 'http', 'headers': [(b'origin', b'http://localhost:3000')]})
        print(_main.settings.ALLOWED_ORIGINS)
        response = await _main.app.exception_handlers[500](request, None)
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://localhost:3000')


class CustomHttpExceptionHandlerTest(unittest.TestCase):
    def create_dummy_client(self, allow_origins):
        dummy_app = FastAPI()
        dummy_app.add_middleware(
            CORSMiddleware,
            allow_origins=allow_origins,
            allow_credentials=True,
            allow_methods=['*'],
            allow_headers=['*'],
            expose_headers=['x-warning'],
        )

        dummy_app.exception_handlers[500] = lambda a, b: custom_http_exception_handler(a, b, dummy_app, allow_origins)

        # Define a route to trigger the exception handler
        @dummy_app.get('/test-trigger-500')
        async def trigger_500(request: Request):
            try:
                raise Exception('Simulated internal server error')
            except Exception as exc:
                # Manually call the custom exception handler
                return await dummy_app.exception_handlers[500](request, exc)

        return TestClient(dummy_app)

    def test_internal_server_error_with_cors_headers(self):
        # Simulate a request from a specific origin
        headers = {'Origin': 'http://example.com'}

        # Origin is not allowed, access-control-allow-origin is not set
        dummy_client = self.create_dummy_client(['http://localhost:3000', 'http://localhost:5173'])
        response = dummy_client.get('/test-trigger-500', headers=headers)
        self.assertEqual(response.status_code, 500)
        self.assertEqual(response.json(), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertNotIn('access-control-allow-origin', response.headers)

        # All origins allowed, it's set to *
        dummy_client = self.create_dummy_client(['*'])
        response = dummy_client.get('/test-trigger-500', headers=headers)
        self.assertEqual(response.status_code, 500)
        self.assertEqual(response.json(), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], '*')

        headers = {'Origin': 'http://localhost:3000'}

        # Origin allowed, return that origin
        dummy_client = self.create_dummy_client(['http://localhost:3000', 'http://localhost:5173'])
        response = dummy_client.get('/test-trigger-500', headers=headers)
        self.assertEqual(response.status_code, 500)
        self.assertEqual(response.json(), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://localhost:3000')

        headers = {'Origin': 'http://localhost:3000'}

        # With cookies, origin returned even it all allowed
        dummy_client = self.create_dummy_client(['*'])
        response = dummy_client.get('/test-trigger-500', headers=headers, cookies={'session': 'abc123'})
        self.assertEqual(response.status_code, 500)
        self.assertEqual(response.json(), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://localhost:3000')
